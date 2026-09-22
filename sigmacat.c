/* sigmacat - nhlt dmic detector + sigma-delta model
 * author: cocomelonc
 *
 * reads an acpi nhlt (non-hda link table), reports every audio endpoint and
 * flags the ones whose linktype is pdm (== a digital microphone, dmic).
 * a dmic is physically a 1-bit sigma-delta modulator feeding a pdm stream,
 * so for each dmic we also run the matching math: analytical peak snr plus a
 * time-domain modulator+decimator simulation that confirms it.
 *
 * deps: libc + libm only.
 * build: gcc -O2 -Wall -Wextra -std=c11 -o sigmacat sigmacat.c -lm
 *
 * usage:
 *   sudo ./sigmacat               # read /sys/firmware/acpi/tables/NHLT
 *   ./sigmacat --raw dump.bin     # read a raw binary nhlt blob
 *   ./sigmacat --hex /tmp/nhlt.hex# read an `xxd` hexdump
 *   ./sigmacat --demo             # no hardware: model a 48k/16-bit dmic
 *   ./sigmacat --viz              # live spectrum of the synthetic model
 *   ./sigmacat --wav rec.wav      # live spectrum + descriptors of real audio
 *   ./sigmacat --scan dir/        # scan a dataset of .wav files
 *   ./sigmacat ... --pdm-clock 3072000
 */
#define _POSIX_C_SOURCE 200809L   /* nanosleep */
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <dirent.h>

#define NHLT_DEFAULT_PATH "/sys/firmware/acpi/tables/NHLT"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- nhlt link types (endpoint linktype field) ------------------------- */
enum { NHLT_HDA = 0, NHLT_DSP = 1, NHLT_PDM = 2, NHLT_SSP = 3,
     NHLT_SLIMBUS = 4, NHLT_SDW = 5 };

static const char *link_name(uint8_t lt) {
  switch (lt) {
  case NHLT_HDA:   return "hd-audio";
  case NHLT_DSP:   return "dsp";
  case NHLT_PDM:   return "pdm (dmic)";
  case NHLT_SSP:   return "ssp (i2s/tdm)";
  case NHLT_SLIMBUS: return "slimbus";
  case NHLT_SDW:   return "soundwire";
  default:       return "unknown";
  }
}

static const char *dir_name(uint8_t d) {
  return d == 0 ? "render" : d == 1 ? "capture" : d == 2 ? "bidir" : "?";
}

/* ---- packed acpi/nhlt layout ------------------------------------------- */
/* standard 36-byte acpi table header; length spans the whole table. */
struct acpi_hdr {
  char   sig[4];
  uint32_t length;
  uint8_t  revision;
  uint8_t  checksum;
  char   oem_id[6];
  char   oem_table_id[8];
  uint32_t oem_revision;
  char   asl_id[4];
  uint32_t asl_revision;
} __attribute__((packed));

/* fixed part of an endpoint descriptor; `length` spans the whole descriptor
 * (fixed part + specific config + formats), so it is how you reach the next. */
struct nhlt_ep {
  uint32_t length;
  uint8_t  link_type;
  uint8_t  instance_id;
  uint16_t vendor_id;
  uint16_t device_id;
  uint16_t revision_id;
  uint32_t subsystem_id;
  uint8_t  device_type;
  uint8_t  direction;
  uint8_t  virtual_bus_id;
} __attribute__((packed));

/* classic 18-byte waveformatex prefix; cb_size tells the extension length. */
struct wav_fmt {
  uint16_t fmt_tag;
  uint16_t channels;
  uint32_t samples_per_sec;
  uint32_t avg_bytes_per_sec;
  uint16_t block_align;
  uint16_t bits_per_sample;
  uint16_t cb_size;
} __attribute__((packed));

/* alignment-safe little-endian reads (we target x86). */
static uint32_t rd32(const uint8_t *p) { uint32_t v; memcpy(&v, p, 4); return v; }

/* =======================================================================
 * sigma-delta / pdm math
 * =====================================================================*/

/* closed-form peak snr of an ideal l-th order, n-bit sigma-delta modulator.
 * the quantization noise is shaped by ntf(z) = (1 - z^-1)^l; integrating
 * |ntf|^2 across the in-band region gives
 *   snr = 6.02*n + 1.76 + 10*log10((2l+1)/pi^(2l)) + 10*(2l+1)*log10(osr).
 */
static double sd_theory_snr(double osr, int order, int nbits) {
  double l = (double)order;
  double shape = 10.0 * log10((2.0 * l + 1.0) / pow(M_PI, 2.0 * l));
  double gain  = 10.0 * (2.0 * l + 1.0) * log10(osr);
  return 6.02 * nbits + 1.76 + shape + gain;
}

/* input sine amplitude relative to the +-1 quantizer full scale. kept below
 * full scale so the 1-bit 2nd-order loop stays stable; the resulting input
 * backoff (in db) is applied to the theory so the two are compared fairly. */
#define SD_INPUT_AMP 0.5

/* iterative radix-2 cooley-tukey fft, in place; n must be a power of two.
 * re[]/im[] hold the complex signal on entry and its spectrum on return. */
static void fft(double *re, double *im, int n) {
  /* bit-reversal permutation */
  for (int i = 1, j = 0; i < n; i++) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      double tr = re[i]; re[i] = re[j]; re[j] = tr;
      double ti = im[i]; im[i] = im[j]; im[j] = ti;
    }
  }
  /* butterflies */
  for (int len = 2; len <= n; len <<= 1) {
    double ang = -2.0 * M_PI / len, wr = cos(ang), wi = sin(ang);
    for (int i = 0; i < n; i += len) {
      double cr = 1.0, ci = 0.0;
      for (int k = 0; k < len / 2; k++) {
        int a = i + k, b = a + len / 2;
        double tr = re[b] * cr - im[b] * ci;
        double ti = re[b] * ci + im[b] * cr;
        re[b] = re[a] - tr; im[b] = im[a] - ti;
        re[a] += tr;        im[a] += ti;
        double ncr = cr * wr - ci * wi;
        ci = cr * wi + ci * wr; cr = ncr;
      }
    }
  }
}

/* simulate a 2nd-order 1-bit modulator on a pure tone, decimate with a
 * sinc^3 (cic) filter, then estimate in-band snr with welch's method:
 * 50%-overlapping hann-windowed segments whose periodograms are averaged.
 * the averaging kills the variance that made a single rectangular dft read
 * ~25 db low. returns measured snr in db; writes tone freq and clamped osr. */
static double sd_simulate(int osr, double f_pcm, double *tone_hz_out, int *osr_out) {
  if (osr < 8)  osr = 8;         /* keep the run small and stable */
  if (osr > 64) osr = 64;
  *osr_out = osr;

  const int seg = 4096;          /* fft / welch segment length (power of 2) */
  const int n_seg = 32;          /* number of 50%-overlap segments to average */
  const int tone_k = 127;        /* signal bin: on-grid so no scalloping loss */
  const int guard = 3;           /* bins each side of the tone (hann mainlobe) */
  const int warm = seg;          /* discard modulator/cic startup transient */

  int n_pcm = warm + seg + (n_seg - 1) * (seg / 2);
  long n_pdm = (long)n_pcm * osr;
  double f_pdm = f_pcm * osr;
  double tone_hz = (double)tone_k / seg * f_pcm;
  *tone_hz_out = tone_hz;

  double *pdm = malloc((size_t)n_pdm * sizeof(double));
  double *a   = malloc((size_t)n_pdm * sizeof(double));
  double *b   = malloc((size_t)n_pdm * sizeof(double));
  double *pcm = malloc((size_t)n_pcm * sizeof(double));
  double *re  = malloc((size_t)seg   * sizeof(double));
  double *im  = malloc((size_t)seg   * sizeof(double));
  double *psd = calloc((size_t)seg / 2, sizeof(double));
  double *win = malloc((size_t)seg   * sizeof(double));
  if (!pdm || !a || !b || !pcm || !re || !im || !psd || !win) {
    free(pdm); free(a); free(b); free(pcm); free(re); free(im); free(psd); free(win);
    return 0.0;
  }

  /* 2nd-order modulator: two integrators, 1-bit quantizer in the loop. */
  double i1 = 0.0, i2 = 0.0, y = 0.0;
  for (long n = 0; n < n_pdm; n++) {
    double x = SD_INPUT_AMP * sin(2.0 * M_PI * tone_hz * n / f_pdm);
    i1 += x - y;
    i2 += i1 - y;
    y = (i2 >= 0.0) ? 1.0 : -1.0;
    pdm[n] = y;
  }

  /* sinc^3 decimator: length-osr moving average cascaded 3x, then pick. */
  memcpy(a, pdm, (size_t)n_pdm * sizeof(double));
  for (int s = 0; s < 3; s++) {
    double acc = 0.0;
    for (long i = 0; i < n_pdm; i++) {
      acc += a[i];
      if (i >= osr) acc -= a[i - osr];
      b[i] = acc / osr;
    }
    double *t = a; a = b; b = t;
  }
  for (int i = 0; i < n_pcm; i++) pcm[i] = a[(long)(i + 1) * osr - 1];

  /* hann window (periodic form matches welch overlap-add analysis). */
  for (int i = 0; i < seg; i++)
    win[i] = 0.5 - 0.5 * cos(2.0 * M_PI * i / seg);

  /* accumulate windowed periodograms over overlapping segments. */
  for (int s = 0; s < n_seg; s++) {
    int start = warm + s * (seg / 2);
    for (int i = 0; i < seg; i++) { re[i] = pcm[start + i] * win[i]; im[i] = 0.0; }
    fft(re, im, seg);
    for (int k = 0; k < seg / 2; k++) psd[k] += re[k] * re[k] + im[k] * im[k];
  }

  /* signal = tone bins +- guard; noise = the rest of the in-band spectrum
   * (skip dc bin 0, which the window leaves any residual mean in). */
  double sig = 0.0, noise = 0.0;
  for (int k = 1; k < seg / 2; k++) {
    if (k >= tone_k - guard && k <= tone_k + guard) sig += psd[k];
    else noise += psd[k];
  }

  free(pdm); free(a); free(b); free(pcm); free(re); free(im); free(psd); free(win);
  return noise > 0.0 ? 10.0 * log10(sig / noise) : 0.0;
}

/* print the full model for one dmic format. */
static void model_dmic(uint32_t f_pcm, uint16_t bits, double f_pdm) {
  double osr = f_pdm / f_pcm;
  int order = 2;              /* typical low-order dmic loop */

  /* full-scale peak snr, then the same theory scaled to the input level the
   * simulation actually drives (the quantization noise floor is set by the
   * modulator/osr, so snr tracks input level db-for-db). */
  double backoff = 20.0 * log10(SD_INPUT_AMP);
  double snr_fs  = sd_theory_snr(osr, order, 1);
  double snr_adj = snr_fs + backoff;

  double tone_hz; int osr_sim;
  double snr_sim = sd_simulate((int)lround(osr), f_pcm, &tone_hz, &osr_sim);

  printf("  sigma-delta model (order l=%d, 1-bit quantizer):\n", order);
  printf("    pcm rate            : %u hz  (%u-bit)\n", f_pcm, bits);
  printf("    pdm clock           : %.0f hz\n", f_pdm);
  printf("    osr                 : %.0f  (f_pdm / f_pcm)\n", osr);
  printf("    theory (full-scale) : %6.1f db   -> enob %.1f bits\n",
         snr_fs, (snr_fs - 1.76) / 6.02);
  printf("    theory (@ %.1f dbfs): %6.1f db   -> enob %.1f bits\n",
         backoff, snr_adj, (snr_adj - 1.76) / 6.02);
  printf("    sim (hann+welch)    : %6.1f db   -> enob %.1f bits (osr=%d, tone %.0f hz)\n",
         snr_sim, (snr_sim - 1.76) / 6.02, osr_sim, tone_hz);
}

/* =======================================================================
 * live terminal visualization
 * =====================================================================*/

/* ansi escapes (plain vt100, no ncurses dependency). */
#define ESC "\x1b["
static void sleep_ms(int ms) {
  struct timespec ts = { ms / 1000, (long)(ms % 1000) * 1000000L };
  nanosleep(&ts, NULL);
}

/* shared bar-spectrum renderer: draw a peak-normalized psd as a terminal
 * spectrum with a db y-axis and a 0..fs/2 frequency x-axis. the column nearest
 * hi_hz is highlighted (tone / dominant peak). caller prints header + footer. */
static void draw_bars(const double *psd, int half, double f_pcm, double hi_hz) {
  const int C = 72, H = 18;
  const double top = 0.0, bot = -108.0, dbrow = (top - bot) / H;
  int barrow[128], hi_col = (int)((hi_hz / (f_pcm / 2.0)) * C);

  double peak = 1e-300;
  for (int k = 1; k < half; k++) if (psd[k] > peak) peak = psd[k];

  for (int c = 0; c < C; c++) {
    int k0 = 1 + c * (half - 1) / C, k1 = 1 + (c + 1) * (half - 1) / C;
    double m = 0.0;
    for (int k = k0; k < k1 && k < half; k++) if (psd[k] > m) m = psd[k];
    double db = 10.0 * log10((m + 1e-300) / peak);
    int r = (int)((top - db) / dbrow + 0.5);
    barrow[c] = r < 0 ? 0 : (r > H ? H : r);
  }

  for (int r = 0; r < H; r++) {
    printf("  %4.0f " ESC "2m|" ESC "0m", top - r * dbrow);
    for (int c = 0; c < C; c++) {
      if (r >= barrow[c]) {
        int hot = (c == hi_col);
        printf("%s%c" ESC "0m", hot ? ESC "36;1m" : ESC "32m", hot ? '#' : '|');
      } else putchar(' ');
    }
    putchar('\n');
  }
  printf("       " ESC "2m+"); for (int c = 0; c < C; c++) putchar('-'); printf(ESC "0m\n");
  printf("       0 hz");
  for (int c = 0; c < C - 12; c++) putchar(' ');
  printf("%.0f hz\n", f_pcm / 2.0);
}

/* animate the welch psd building up segment by segment: a terminal spectrum
 * analyzer of the decimated dmic output, with the theoretical noise-shaping
 * asymptote (in-band quantization noise ~ f^(2l)) overlaid and a running snr
 * estimate that converges as more segments are averaged in. */
static void viz_run(int osr, double f_pcm) {
  if (osr < 8)  osr = 8;
  if (osr > 64) osr = 64;

  const int seg = 2048, half = seg / 2, tone_k = 64, guard = 3, order = 2;
  const int n_seg = 60, warm = seg;

  int n_pcm = warm + seg + (n_seg - 1) * (seg / 2);
  long n_pdm = (long)n_pcm * osr;
  double f_pdm = f_pcm * osr;
  double tone_hz = (double)tone_k / seg * f_pcm;

  double *pdm = malloc((size_t)n_pdm * sizeof(double));
  double *a   = malloc((size_t)n_pdm * sizeof(double));
  double *b   = malloc((size_t)n_pdm * sizeof(double));
  double *pcm = malloc((size_t)n_pcm * sizeof(double));
  double *re  = malloc((size_t)seg   * sizeof(double));
  double *im  = malloc((size_t)seg   * sizeof(double));
  double *psd = calloc((size_t)half, sizeof(double));
  double *win = malloc((size_t)seg   * sizeof(double));
  if (!pdm || !a || !b || !pcm || !re || !im || !psd || !win) {
    free(pdm); free(a); free(b); free(pcm); free(re); free(im); free(psd); free(win);
    fprintf(stderr, "viz: out of memory\n"); return;
  }

  /* modulator + sinc^3 decimator (same chain as the headless sim). */
  double i1 = 0.0, i2 = 0.0, y = 0.0;
  for (long n = 0; n < n_pdm; n++) {
    double x = SD_INPUT_AMP * sin(2.0 * M_PI * tone_hz * n / f_pdm);
    i1 += x - y; i2 += i1 - y; y = (i2 >= 0.0) ? 1.0 : -1.0; pdm[n] = y;
  }
  memcpy(a, pdm, (size_t)n_pdm * sizeof(double));
  for (int s = 0; s < 3; s++) {
    double acc = 0.0;
    for (long i = 0; i < n_pdm; i++) { acc += a[i]; if (i >= osr) acc -= a[i - osr]; b[i] = acc / osr; }
    double *t = a; a = b; b = t;
  }
  for (int i = 0; i < n_pcm; i++) pcm[i] = a[(long)(i + 1) * osr - 1];
  for (int i = 0; i < seg; i++) win[i] = 0.5 - 0.5 * cos(2.0 * M_PI * i / seg);

  double snr_fs  = sd_theory_snr(osr, order, 1);
  double snr_adj = snr_fs + 20.0 * log10(SD_INPUT_AMP);   /* at sim input level */

  printf(ESC "?25l" ESC "2J");              /* hide cursor, clear */
  for (int s = 0; s < n_seg; s++) {
    /* fold one more hann-windowed segment into the averaged periodogram. */
    int start = warm + s * (seg / 2);
    for (int i = 0; i < seg; i++) { re[i] = pcm[start + i] * win[i]; im[i] = 0.0; }
    fft(re, im, seg);
    for (int k = 0; k < half; k++) psd[k] += re[k] * re[k] + im[k] * im[k];

    /* running in-band snr from the accumulated psd. */
    double sig = 0.0, noise = 0.0;
    for (int k = 1; k < half; k++) {
      if (k >= tone_k - guard && k <= tone_k + guard) sig += psd[k]; else noise += psd[k];
    }
    double snr = noise > 0.0 ? 10.0 * log10(sig / noise) : 0.0;

    /* draw frame: header + shared bar renderer + footer. */
    printf(ESC "H");
    printf("  sigmacat :: live pdm/sigma-delta spectrum (welch averaging)\n");
    printf("  osr=%d  f_pcm=%.0f hz  f_pdm=%.0f hz  tone=%.0f hz  segment %2d/%d\n\n",
           osr, f_pcm, f_pdm, tone_hz, s + 1, n_seg);
    draw_bars(psd, half, f_pcm, tone_hz);
    putchar('\n');

    printf("  " ESC "32m|" ESC "0m = measured psd   " ESC "36;1m#" ESC "0m = test tone"
           "   (rising floor toward the edge = sigma-delta noise shaping)\n");
    printf("  in-band snr : " ESC "1m%6.1f db" ESC "0m (enob %.1f)   converging to "
           "theory %.1f db (@ %.0f dbfs)  /  %.1f db (full scale)\n",
           snr, (snr - 1.76) / 6.02, snr_adj, 20.0 * log10(SD_INPUT_AMP), snr_fs);
    fflush(stdout);
    sleep_ms(45);
  }
  printf(ESC "?25h");                        /* show cursor */

  free(pdm); free(a); free(b); free(pcm); free(re); free(im); free(psd); free(win);
}

/* =======================================================================
 * real-data spectrum: wav reader, live viewer, dataset scan
 * =====================================================================*/

static uint16_t rd16le(const uint8_t *p) { return (uint16_t)(p[0] | (p[1] << 8)); }
static int cmp_d(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}
static int is_wav(const char *nm) {
  size_t l = strlen(nm);
  if (l < 4 || nm[l - 4] != '.') return 0;
  return tolower(nm[l-3]) == 'w' && tolower(nm[l-2]) == 'a' && tolower(nm[l-1]) == 'v';
}

/* read a pcm wav (16-/32-bit int or 32-bit float; any channel count) into a
 * malloc'd mono double buffer normalized to ~[-1,1]. returns per-channel frame
 * count; writes sample rate and channel count out. */
static double *wav_read(const char *path, long *n_out, int *sr_out, int *ch_out) {
  FILE *f = fopen(path, "rb");
  if (!f) { perror(path); return NULL; }

  uint8_t hdr[12];
  if (fread(hdr, 1, 12, f) != 12 || memcmp(hdr, "RIFF", 4) || memcmp(hdr + 8, "WAVE", 4)) {
    fprintf(stderr, "%s: not a RIFF/WAVE file\n", path); fclose(f); return NULL;
  }

  int fmt = 0, ch = 0, sr = 0, bits = 0;
  double *out = NULL; long nframes = 0;
  uint8_t ch_hdr[8];

  while (fread(ch_hdr, 1, 8, f) == 8) {
    uint32_t sz = rd32(ch_hdr + 4);
    if (!memcmp(ch_hdr, "fmt ", 4)) {
      uint8_t fb[40];
      uint32_t rd = sz > sizeof fb ? (uint32_t)sizeof fb : sz;
      if (fread(fb, 1, rd, f) != rd) break;
      fmt  = rd16le(fb);
      ch   = rd16le(fb + 2);
      sr   = (int)rd32(fb + 4);
      bits = rd16le(fb + 14);
      if (fmt == 0xFFFE && rd >= 26) fmt = rd16le(fb + 24);   /* extensible subformat */
      if (sz > rd) fseek(f, (long)(sz - rd), SEEK_CUR);
    } else if (!memcmp(ch_hdr, "data", 4)) {
      if (ch <= 0 || sr <= 0 || (bits != 16 && bits != 32)) {
        fprintf(stderr, "%s: unsupported format (fmt=%d ch=%d bits=%d)\n", path, fmt, ch, bits);
        break;
      }
      int bytes = bits / 8;
      nframes = (long)sz / bytes / ch;
      uint8_t *raw = malloc(sz);
      out = malloc((size_t)nframes * sizeof(double));
      if (!raw || !out || fread(raw, 1, sz, f) != sz) { free(raw); free(out); out = NULL; break; }
      for (long i = 0; i < nframes; i++) {
        double acc = 0.0;                     /* downmix channels to mono */
        for (int c = 0; c < ch; c++) {
          uint8_t *p = raw + ((long)i * ch + c) * bytes;
          double v;
          if (bits == 16)      { int16_t s; memcpy(&s, p, 2); v = s / 32768.0; }
          else if (fmt == 3)   { float   s; memcpy(&s, p, 4); v = s; }
          else                 { int32_t s; memcpy(&s, p, 4); v = s / 2147483648.0; }
          acc += v;
        }
        out[i] = acc / ch;
      }
      free(raw);
      break;                                  /* first data chunk is enough */
    } else {
      fseek(f, (long)(sz + (sz & 1)), SEEK_CUR);  /* skip unknown chunk (+pad) */
    }
  }
  fclose(f);
  if (!out) return NULL;
  *n_out = nframes; *sr_out = sr; *ch_out = ch;
  return out;
}

/* averaged welch psd over a real mono signal (hann, 50% overlap). */
static int welch_psd(const double *x, long n, int seg, double *psd) {
  int half = seg / 2, nseg = 0;
  double *re = malloc((size_t)seg * sizeof(double));
  double *im = malloc((size_t)seg * sizeof(double));
  double *win = malloc((size_t)seg * sizeof(double));
  if (!re || !im || !win) { free(re); free(im); free(win); return 0; }
  for (int i = 0; i < seg; i++) win[i] = 0.5 - 0.5 * cos(2.0 * M_PI * i / seg);
  for (int k = 0; k < half; k++) psd[k] = 0.0;
  for (long start = 0; start + seg <= n; start += seg / 2) {
    double mean = 0.0;                         /* per-segment dc removal */
    for (int i = 0; i < seg; i++) mean += x[start + i];
    mean /= seg;
    for (int i = 0; i < seg; i++) { re[i] = (x[start + i] - mean) * win[i]; im[i] = 0.0; }
    fft(re, im, seg);
    for (int k = 0; k < half; k++) psd[k] += re[k] * re[k] + im[k] * im[k];
    nseg++;
  }
  if (nseg) for (int k = 0; k < half; k++) psd[k] /= nseg;
  free(re); free(im); free(win);
  return nseg;
}

/* descriptors of a captured signal. */
struct sig_stats {
  double rms_dbfs, crest_db, dom_hz, centroid_hz, floor_db, flatness_db;
};

/* time- and frequency-domain descriptors from the signal and its psd. */
static struct sig_stats compute_stats(const double *x, long n, int sr,
                                      const double *psd, int half, int seg) {
  struct sig_stats st; memset(&st, 0, sizeof st);

  double sumsq = 0.0, peak = 0.0;
  for (long i = 0; i < n; i++) { double a = fabs(x[i]); if (a > peak) peak = a; sumsq += x[i] * x[i]; }
  double rms = sqrt(sumsq / n);
  st.rms_dbfs = rms > 0 ? 20.0 * log10(rms) : -120.0;
  st.crest_db = (rms > 0 && peak > 0) ? 20.0 * log10(peak / rms) : 0.0;

  /* analyze [~50 hz, fs/2] so residual dc/rumble does not skew descriptors. */
  int kmin = (int)(50.0 * seg / sr); if (kmin < 1) kmin = 1;
  int m = half - kmin;
  double pk = 1e-300, sump = 0.0, sumfp = 0.0, sumlog = 0.0, suma = 0.0;
  int dom = kmin;
  for (int k = kmin; k < half; k++) {
    if (psd[k] > pk) { pk = psd[k]; dom = k; }
    double fk = (double)k * sr / seg;
    sump += psd[k]; sumfp += fk * psd[k];
    sumlog += log(psd[k] + 1e-300); suma += psd[k];
  }
  st.dom_hz = (double)dom * sr / seg;
  st.centroid_hz = sump > 0 ? sumfp / sump : 0.0;
  st.flatness_db = (m > 0 && suma > 0) ? 10.0 * log10(exp(sumlog / m) / (suma / m)) : 0.0;

  double *tmp = malloc((size_t)m * sizeof(double));
  if (tmp) {
    for (int k = kmin; k < half; k++) tmp[k - kmin] = psd[k];
    qsort(tmp, m, sizeof(double), cmp_d);
    st.floor_db = 10.0 * log10((tmp[m / 2] + 1e-300) / pk);   /* median below peak */
    free(tmp);
  }
  return st;
}

/* --wav: live welch spectrum of a real recording + signal descriptors. */
static void wav_view(const char *path) {
  long n; int sr, ch;
  double *x = wav_read(path, &n, &sr, &ch);
  if (!x) return;

  const int seg = 2048, half = seg / 2;
  if (n < seg) { fprintf(stderr, "%s: too short (%ld samples)\n", path, n); free(x); return; }

  double *re  = malloc((size_t)seg * sizeof(double));
  double *im  = malloc((size_t)seg * sizeof(double));
  double *win = malloc((size_t)seg * sizeof(double));
  double *psd = calloc((size_t)half, sizeof(double));
  if (!re || !im || !win || !psd) { free(x); free(re); free(im); free(win); free(psd); return; }
  for (int i = 0; i < seg; i++) win[i] = 0.5 - 0.5 * cos(2.0 * M_PI * i / seg);

  const char *base = strrchr(path, '/'); base = base ? base + 1 : path;
  long total = (n - seg) / (seg / 2) + 1, done = 0;
  int every = (int)(total / 60); if (every < 1) every = 1;

  printf(ESC "?25l" ESC "2J");
  for (long start = 0; start + seg <= n; start += seg / 2) {
    double mean = 0.0;                         /* per-segment dc removal */
    for (int i = 0; i < seg; i++) mean += x[start + i];
    mean /= seg;
    for (int i = 0; i < seg; i++) { re[i] = (x[start + i] - mean) * win[i]; im[i] = 0.0; }
    fft(re, im, seg);
    for (int k = 0; k < half; k++) psd[k] += re[k] * re[k] + im[k] * im[k];
    done++;

    if (done % every == 0 || start + seg + seg / 2 > n) {
      struct sig_stats st = compute_stats(x, n, sr, psd, half, seg);
      printf(ESC "H");
      printf("  sigmacat :: real-data spectrum  " ESC "1m%s" ESC "0m\n", base);
      printf("  %d hz  %d ch  %.2f s  |  welch %ld/%ld hann segments\n\n",
             sr, ch, (double)n / sr, done, total);
      draw_bars(psd, half, (double)sr, st.dom_hz);
      putchar('\n');
      printf("  " ESC "36;1m#" ESC "0m dominant " ESC "1m%.0f hz" ESC "0m   "
             "centroid %.0f hz   rms %.1f dbfs   crest %.1f db\n",
             st.dom_hz, st.centroid_hz, st.rms_dbfs, st.crest_db);
      printf("  noise floor %.1f db below peak   spectral flatness %.1f db (%s)\n",
             st.floor_db, st.flatness_db, st.flatness_db > -10.0 ? "noise-like" : "tonal");
      fflush(stdout);
      sleep_ms(45);
    }
  }
  printf(ESC "?25h");
  free(x); free(re); free(im); free(win); free(psd);
}

/* --scan: run the analysis over every .wav in a directory (a dataset) and
 * print one descriptor row per file plus aggregate means. */
static void scan_dataset(const char *dir) {
  DIR *d = opendir(dir);
  if (!d) { perror(dir); return; }

  const int seg = 2048, half = seg / 2;
  double *psd = malloc((size_t)half * sizeof(double));
  if (!psd) { closedir(d); return; }

  printf("dataset scan: %s\n\n", dir);
  printf("%-24s %6s %3s %7s %8s %8s %8s %8s\n",
         "file", "sr", "ch", "dur_s", "dom_hz", "cent_hz", "rms_dbfs", "flat_db");
  printf("%-24s %6s %3s %7s %8s %8s %8s %8s\n",
         "------------------------", "------", "---", "-------",
         "--------", "--------", "--------", "--------");

  struct dirent *e; int files = 0; double sdom = 0, scent = 0, srms = 0;
  while ((e = readdir(d))) {
    if (!is_wav(e->d_name)) continue;
    char full[4096];
    snprintf(full, sizeof full, "%s/%s", dir, e->d_name);
    long n; int sr, ch;
    double *x = wav_read(full, &n, &sr, &ch);
    if (!x) continue;
    if (n < seg) { free(x); continue; }
    welch_psd(x, n, seg, psd);
    struct sig_stats st = compute_stats(x, n, sr, psd, half, seg);
    printf("%-24.24s %6d %3d %7.2f %8.0f %8.0f %8.1f %8.1f\n",
           e->d_name, sr, ch, (double)n / sr,
           st.dom_hz, st.centroid_hz, st.rms_dbfs, st.flatness_db);
    files++; sdom += st.dom_hz; scent += st.centroid_hz; srms += st.rms_dbfs;
    free(x);
  }
  free(psd); closedir(d);

  if (files)
    printf("\n%d file%s  |  mean: dominant %.0f hz  centroid %.0f hz  rms %.1f dbfs\n",
           files, files == 1 ? "" : "s", sdom / files, scent / files, srms / files);
  else
    printf("\nno .wav files found in %s\n", dir);
}

/* =======================================================================
 * nhlt parsing
 * =====================================================================*/

/* walk one endpoint's formats block and model each pdm format. */
static void handle_formats(const uint8_t *p, const uint8_t *end, double f_pdm) {
  if (p >= end) return;
  uint8_t count = *p++;

  for (uint8_t i = 0; i < count; i++) {
    if (p + sizeof(struct wav_fmt) > end) break;
    struct wav_fmt wf;
    memcpy(&wf, p, sizeof wf);

    printf("  format[%u]  : %u hz, %u ch, %u-bit\n",
         i, wf.samples_per_sec, wf.channels, wf.bits_per_sample);
    model_dmic(wf.samples_per_sec, wf.bits_per_sample, f_pdm);

    /* advance: waveformatex + extension + nhlt_specific_cfg(u32 size+data) */
    const uint8_t *q = p + sizeof(struct wav_fmt) + wf.cb_size;
    if (q + 4 > end) break;
    q += 4 + rd32(q);
    if (q <= p || q > end) break;
    p = q;
  }
}

static int parse_nhlt(const uint8_t *buf, size_t len, double f_pdm) {
  if (len < sizeof(struct acpi_hdr)) {
    fprintf(stderr, "nhlt: too small for acpi header (%zu bytes)\n", len);
    return 1;
  }
  const struct acpi_hdr *h = (const struct acpi_hdr *)buf;
  if (memcmp(h->sig, "NHLT", 4) != 0) {
    fprintf(stderr, "nhlt: bad signature '%.4s' (expected NHLT)\n", h->sig);
    return 1;
  }
  if (h->length > len) {
    fprintf(stderr, "nhlt: header length %u > buffer %zu\n", h->length, len);
    return 1;
  }

  printf("nhlt table   : %u bytes, oem '%.6s', rev %u\n",
       h->length, h->oem_id, h->revision);

  const uint8_t *p = buf + sizeof(struct acpi_hdr);
  const uint8_t *end = buf + h->length;
  if (p >= end) { fprintf(stderr, "nhlt: no endpoint count\n"); return 1; }

  uint8_t ep_count = *p++;
  printf("endpoints  : %u\n", ep_count);

  int dmic_total = 0;
  for (int i = 0; i < ep_count; i++) {
    if (p + sizeof(struct nhlt_ep) > end) {
      fprintf(stderr, "nhlt: truncated endpoint %d\n", i); return 1;
    }
    struct nhlt_ep ep;
    memcpy(&ep, p, sizeof ep);
    if (ep.length < sizeof(struct nhlt_ep) || p + ep.length > end) {
      fprintf(stderr, "nhlt: endpoint %d bad length %u\n", i, ep.length);
      return 1;
    }

    int is_dmic = (ep.link_type == NHLT_PDM);
    printf("\n[ep %d] linktype=%u %-13s dir=%s vid=%04x did=%04x%s\n",
         i, ep.link_type, link_name(ep.link_type),
         dir_name(ep.direction), ep.vendor_id, ep.device_id,
         is_dmic ? "   <== DMIC" : "");

    if (is_dmic) {
      dmic_total++;
      const uint8_t *body = p + sizeof(struct nhlt_ep);
      const uint8_t *ep_end = p + ep.length;
      if (body + 4 <= ep_end) {
        const uint8_t *fmts = body + 4 + rd32(body); /* skip specific cfg */
        if (fmts < ep_end) handle_formats(fmts, ep_end, f_pdm);
      }
    }
    p += ep.length;
  }

  printf("\nverdict    : %s (%d dmic endpoint%s declared in nhlt)\n",
       dmic_total ? "DMIC PRESENT" : "no dmic",
       dmic_total, dmic_total == 1 ? "" : "s");
  return 0;
}

/* =======================================================================
 * loaders + cli
 * =====================================================================*/

/* slurp a whole file into a malloc'd buffer. */
static uint8_t *read_file(const char *path, size_t *out_len) {
  FILE *f = fopen(path, "rb");
  if (!f) { perror(path); return NULL; }
  fseek(f, 0, SEEK_END);
  long sz = ftell(f);
  fseek(f, 0, SEEK_SET);
  if (sz <= 0) { fclose(f); fprintf(stderr, "%s: empty\n", path); return NULL; }

  uint8_t *buf = malloc((size_t)sz);
  size_t got = buf ? fread(buf, 1, (size_t)sz, f) : 0;
  fclose(f);
  if (!buf || got != (size_t)sz) { free(buf); return NULL; }
  *out_len = got;
  return buf;
}

/* decode an `xxd` hexdump: take the two-hex-digit tokens between the offset
 * colon and the ascii gutter (two-or-more spaces). tolerant of plain hex too. */
static uint8_t *read_hexdump(const char *path, size_t *out_len) {
  FILE *f = fopen(path, "rb");
  if (!f) { perror(path); return NULL; }

  size_t cap = 4096, n = 0;
  uint8_t *buf = malloc(cap);
  char line[512];
  if (!buf) { fclose(f); return NULL; }

  while (fgets(line, sizeof line, f)) {
    char *s = strchr(line, ':');       /* skip offset field if present */
    s = s ? s + 1 : line;

    while (*s) {
      /* stop at the ascii gutter: two spaces in a row after some hex */
      if (s[0] == ' ' && s[1] == ' ') break;
      if (isxdigit((unsigned char)s[0]) && isxdigit((unsigned char)s[1])) {
        unsigned v;
        sscanf(s, "%2x", &v);
        if (n == cap) { cap *= 2; buf = realloc(buf, cap); if (!buf) { fclose(f); return NULL; } }
        buf[n++] = (uint8_t)v;
        s += 2;
      } else {
        s++;
      }
    }
  }
  fclose(f);
  *out_len = n;
  return buf;
}

static void usage(const char *argv0) {
  fprintf(stderr,
    "sigmacat - nhlt dmic detector + sigma-delta model (author: cocomelonc)\n"
    "usage:\n"
    "  %s             read %s\n"
    "  %s --raw <file>      raw binary nhlt blob\n"
    "  %s --hex <file>      xxd hexdump\n"
    "  %s --demo          no hardware: model a 48k/16-bit dmic\n"
    "  %s --viz           live spectrum animation of the synthetic model\n"
    "  %s --wav <file>    live spectrum + descriptors of a real recording\n"
    "  %s --scan <dir>    scan a dataset of .wav files, print a metrics table\n"
    "  %s [...] --pdm-clock <hz>  assumed pdm bit clock (default 3072000)\n",
    argv0, NHLT_DEFAULT_PATH, argv0, argv0, argv0, argv0, argv0, argv0, argv0);
}

int main(int argc, char **argv) {
  const char *path = NHLT_DEFAULT_PATH;
  int mode_hex = 0, mode_demo = 0, mode_viz = 0, use_default = 1;
  const char *wav_path = NULL, *scan_dir = NULL;
  double f_pdm = 3072000.0;         /* common dmic clock */

  for (int i = 1; i < argc; i++) {
    if    (!strcmp(argv[i], "--raw") && i + 1 < argc) { path = argv[++i]; use_default = 0; }
    else if (!strcmp(argv[i], "--hex") && i + 1 < argc) { path = argv[++i]; mode_hex = 1; use_default = 0; }
    else if (!strcmp(argv[i], "--demo")) mode_demo = 1;
    else if (!strcmp(argv[i], "--viz")) mode_viz = 1;
    else if (!strcmp(argv[i], "--wav") && i + 1 < argc) wav_path = argv[++i];
    else if (!strcmp(argv[i], "--scan") && i + 1 < argc) scan_dir = argv[++i];
    else if (!strcmp(argv[i], "--pdm-clock") && i + 1 < argc) f_pdm = atof(argv[++i]);
    else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage(argv[0]); return 0; }
    else { fprintf(stderr, "unknown arg: %s\n", argv[i]); usage(argv[0]); return 2; }
  }

  if (wav_path)  { wav_view(wav_path);   return 0; }
  if (scan_dir)  { scan_dataset(scan_dir); return 0; }

  if (mode_viz) {
    viz_run((int)lround(f_pdm / 48000.0), 48000.0);
    return 0;
  }

  if (mode_demo) {
    printf("demo mode  : no nhlt, modeling a synthetic 48000 hz / 16-bit dmic\n\n");
    printf("[ep 0] synthetic %s   <== DMIC\n", link_name(NHLT_PDM));
    printf("  format[0]  : 48000 hz, 2 ch, 16-bit\n");
    model_dmic(48000, 16, f_pdm);
    printf("\nverdict    : DMIC PRESENT (1 synthetic dmic endpoint)\n");
    return 0;
  }

  size_t len = 0;
  uint8_t *buf = mode_hex ? read_hexdump(path, &len) : read_file(path, &len);
  if (!buf) {
    if (use_default)
      fprintf(stderr, "hint: reading %s needs root; try `sudo %s`, or --demo\n",
          NHLT_DEFAULT_PATH, argv[0]);
    return 1;
  }

  int rc = parse_nhlt(buf, len, f_pdm);
  free(buf);
  return rc;
}

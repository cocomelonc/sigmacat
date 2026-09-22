# sigmacat

`author: cocomelonc`

single-file C tool that parses the ACPI *NHLT* (Non-HDA Link Table), reports every audio endpoint, and flags digital microphones - plus the math behind what a DMIC actually is.

## what it does

**detect DMIC.** in NHLT every endpoint carries a `LinkType` byte. a DMIC is a **PDM** link (`LinkType == 2`). the parser follows the NHLT layout by exact offsets - the 36-byte ACPI header, the endpoint-count byte, then each endpoint descriptor whose `Length` field (offset 0) spans the whole descriptor, so you hop endpoint->endpoint by it. it also decodes the advertised audio format(s) (rate / channels / depth).

**model the DMIC.** a DMIC is physically a *1-bit sigma-delta (ΔΣ) modulator* emitting a PDM stream that is decimated to PCM. from the format NHLT declares, sigmacat computes:

- oversampling ratio $\mathrm{OSR} = f_\mathrm{pdm} / f_\mathrm{pcm}$
- closed-form peak SNR of an ideal $L$-th order, $N$-bit modulator:

$$
\mathrm{SNR}_\mathrm{dB} = 6.02\,N + 1.76 + 10\log_{10}\!\left(\frac{2L+1}{\pi^{2L}}\right) + 10\,(2L+1)\log_{10}(\mathrm{OSR})
$$

  which follows from shaping the quantization noise by $\mathrm{NTF}(z) = (1 - z^{-1})^{L}$ and integrating $|\mathrm{NTF}(e^{j\omega})|^2$ over the in-band region $[0, f_\mathrm{pcm}/2]$. the effective resolution is then

$$
\mathrm{ENOB} = \frac{\mathrm{SNR}_\mathrm{dB} - 1.76}{6.02}
$$

- a **live time-domain simulation**: a 2nd-order 1-bit modulator on a test tone $\rightarrow$ sinc³ (CIC) decimator $\rightarrow$ single-bin DFT that measures the real SNR and confirms the analytical number.

how to build:      

```bash
make            # cc -O2 -o sigmacat sigmacat.c -lm  (libc + libm only)
```

run

```bash
sudo ./sigmacat                 # read /sys/firmware/acpi/tables/NHLT (root)
./sigmacat --raw dump.bin       # raw binary nhlt blob
./sigmacat --hex /tmp/nhlt.hex  # an `xxd` hexdump
./sigmacat --demo               # no hardware: model a 48k/16-bit dmic
./sigmacat ... --pdm-clock 3072000
```

original firmware-forensics flow:     

```bash
sudo xxd /sys/firmware/acpi/tables/NHLT > /tmp/nhlt.hex
./sigmacat --hex /tmp/nhlt.hex
```

## notes

the analytical SNR is the *ideal* peak; the built-in simulation runs a deliberately simple 2nd-order loop with a rectangular-window DFT, so its measured SNR sits below theory (window leakage + non-ideal decimation). the gap is expected and is where the deeper DSP pass + live visualization go next.     

LinkType codes: `0` hd-audio, `1` dsp, `2` *pdm/dmic*, `3` ssp/i2s, `4` slimbus, `5` soundwire.

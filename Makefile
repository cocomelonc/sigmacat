# sigmacat - author: cocomelonc
CC     ?= cc
CFLAGS ?= -O2 -Wall -Wextra -std=c11
LDLIBS := -lm

sigmacat: sigmacat.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f sigmacat

.PHONY: clean

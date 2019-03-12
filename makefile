#
# Makefile for kmer_freq C binary
#
# Compiler: gcc 7.3.0
#

HDRS = src/kseq.h src/kmer_freq.h
SRCS = src/kmer_freq.c
CC = gcc
CFLAGS = -lm

kmer_freq: $(SRCS) $(HDRS)
	$(CC) $(SRCS) -o kmer_freq $(CFLAGS)



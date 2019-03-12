#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include "kseq.h"

#define KMER_MIN_SIZE 2
#define KMER_MAX_SIZE 10
#define RC_TERM -11

typedef struct {
	char kmerSize;
	int size;
	int* kmerCountArray;

} KmerCountMap;

void calcKmerCounts(char*, int, KmerCountMap*);
int incKmerCount(char*, int, int, KmerCountMap*);
int getKmerIndex(char*, int, int);
void indexToKmer(int, int, char*);
char getBase(int);
char getCompliment(char);
int initKmerCountMap(int, KmerCountMap*);
void resetKmerCountMap(KmerCountMap*);
void deleteKmerCountMap(KmerCountMap*);
void dumpResults(KmerCountMap*);
void printKmers(KmerCountMap*);
void printFreqs(KmerCountMap*);

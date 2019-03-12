#include "kmer_freq.h"
//macro to define use of file descriptor and read function
KSEQ_INIT(int, read)



/*
 * kmer_freq
 *
 * Usage: kmer_freq <file> <kmer_size>
 * Inputs: a fasta/q file of sequences, the size of kmers to find
 * Outputs: tab format of frequencies of each kmer. Each sequence has frequencies local to that sequence. Outputs multiple sequences.
 *
 * kmer_size must be in the bounds of 2 and 10 inclusively.
 *
 * Complimentary kmers are counted the same, that is for example TTGC is counted as AACG.
 * */

int main(int argc, char** argv) {
	KmerCountMap kmerCountMap;
	FILE* fp;
	kseq_t* seq;
	int rc, fd, kmerSize;

	if (argc != 3) {
		fprintf(stderr, "Usage: %s <file> <kmer_size>\n", argv[0]);
		return 1;
	}

	kmerSize = atoi(argv[2]); //0 if invalid string
	if (kmerSize < 2 || kmerSize > 10) {
		fprintf(stderr, "kmer_size must be between 2 and 10 inclusive");
		return 2;
	}
	
	initKmerCountMap(kmerSize, &kmerCountMap); //initialize kmerCountMap

	fp = fopen(argv[1], "r");
	fd = fileno(fp);

	seq = kseq_init(fd); //initialize seq

	printf("seq_id\t");
	printKmers(&kmerCountMap); //print kmer possibilities
	printf("\n");	
		
	while ((rc = kseq_read(seq)) >= 0) {
		calcKmerCounts(seq->seq.s, kmerSize, &kmerCountMap);
		
		printf("%s\t", seq->name.s);
		printFreqs(&kmerCountMap);
		printf("\n");
	
		resetKmerCountMap(&kmerCountMap);
	}
	
	kseq_destroy(seq); //destroy seq
	fclose(fp);
		
	return 0;
}

/*
 * Calculates kmer counts of given size for the given sequence string and stores results in given map.
 * */
void calcKmerCounts(char* seqString, int kmerSize, KmerCountMap* kmerCountMap) {

	int i = 0, incRC;
	char base;
	while(1) { //breaks when incKmerCount sees null terminator
		incRC = incKmerCount(seqString, i, i+kmerSize, kmerCountMap);
		
		if (incRC == -11) //means null terminator reached in this kmer
			break;

		if (incRC < 0) {//means bad base in the kmer at index incRC
			i += -incRC + 1; // move seq index past invalid base	
		} else {
			i++;
		}
	}	

}

/*
 * Increments the count stored for kmer given in the seq bounded by kmerStart and kmerEnd exclusively
 *
 * Returns 0 if successful
 * */
int incKmerCount(char* seq, int kmerStart, int kmerEnd, KmerCountMap* kmerMap) {
	int kmerIndex;
	assert(kmerEnd - kmerStart == kmerMap->kmerSize); //prevent memory access outside array
	
	if((kmerIndex = getKmerIndex(seq, kmerStart, kmerEnd)) >= 0) {
		kmerMap->kmerCountArray[kmerIndex] += 1;
		return 0;
	} else {
		return kmerIndex;
	}

}

/*
 *	Return the index corresponding to the given sequence of nucleotides
 *	bounded by kmerStart and kmerEnd exclusively.
 *
 *	If a kmer substring begins with G or T, then we count its reverse instead.
 *	This means complimentary kmers are counted the same. 
 *	
 *	Assumes the bounds are actually correct.
 *	
 *	Returns proper index for given Kmer on success
 *	Returns -i if invalid character (not in {A, G, C, T} case insensitive), where i is index into the kmer substring
 *	Returns -11 if null terminator reached
 * */
int getKmerIndex(char* seq, int kmerStart, int kmerEnd) {

	assert(kmerEnd - kmerStart > 0); //a size 0 kmer is meaningless
	 
	char base; // the base (A/G/C/T) at i
	int index = 0; // the index in the array where the kmer should be; what we should return if successful
	int place = 1; // place in the kmer as a base 4 number (1s, 4s, 16s, etc)
	//reverse each base if ends with G or T case-insensititve
	//since start base has highest place, want base to be lower value between compliment 
	int reverse = seq[kmerStart] == 'T' || seq[kmerStart] == 'C' ||
		      seq[kmerStart] == 't' || seq[kmerStart] == 'c';
	
	int i;
	for (i = kmerEnd - 1; i >= kmerStart; i--) {
		base = reverse? getCompliment(seq[i]) : seq[i];

		switch(base) {
			case 'a':
			case 'A': break; // (adding 0*place as A represented by 0 base 4)
		
			case 'g':
			case 'G':
				index += place; // adding 1*place as C represented by 1 base 4)
				break;
			
			case 'c':
			case 'C':
				index += place * 2; //hopefully its clear by now
				break;

			case 't':
			case 'T':
				index += place * 3;
				break;
			case '\0':
				return RC_TERM; //null terminator
			default:
				return -i;	// this is the negated index in the kmer
			       			// substring where we find a bad character 	
		}

		place *= 4; //since we are in base 4 (we can ensure this is optimized by shifting left twice instead of mult)
	
	}

	return index;

}

/*
 * Stores the kmer string representing the given index into kmer.
 *
 * The index is a base 10 representation of the base 4 kmer,
 * where the base 4 digits are A=0, C=1, G=2, T=3 
 * So this is like converting a base 10 number to base 4.
 *
 * Assumes kmer to store string in is array of size kmerSize + 1.
 * */
void indexToKmer(int index, int kmerSize, char* kmer) {
	int value = index;
	int place = (int) pow(4, kmerSize-1); //first value we subtract (highest fours place)
	int i, digit = 0, mod = 1;
	for(i = 0; i < kmerSize; i++) {
		if (value - place >= 0) { //this means our digit is too small of value
			value = value - place;
			digit++;			
			i--; //have to do another iteration of loop
		} else { //right digit value, so reset digit, and put correct base at right position in kmer
			kmer[i] = getBase(digit);
			digit = 0;
			place /= 4; //divide place by 4 to get next place
		}
	}
	kmer[kmerSize] = '\0'; //end with null terminator
}

/*
 * Returns the base associated with the given value in {0, 1, 2, 3}
 * */
char getBase(int value) {
	assert(value >= 0 && value <= 3);
	switch(value) {
		case 0: return 'A';
		case 1: return 'G';
		case 2: return 'C';
		case 3: return 'T';
	}
}

char getCompliment(char base) {
	switch(base) {
		case 'a':
		case 'A': return 'T';
		case 't':
		case 'T': return 'A';
		case 'g':
		case 'G': return 'C';
		case 'c':
		case 'C': return 'G';
		default: return base;
	}
}

/*
 * Generate array on the heap and store addr in kmerCountArray
 * Return the size of the array?
 * */
int initKmerCountMap(int kmerSize, KmerCountMap* kmerCountMap) {
	assert(kmerSize >= 2 && kmerSize <= 10);
	int kmerArraySize = pow(2, (double) 2*kmerSize - 1);
	kmerCountMap->size = kmerArraySize;
	kmerCountMap->kmerSize = kmerSize;
	kmerCountMap->kmerCountArray = (int*) calloc(kmerArraySize, sizeof(int));
	return kmerArraySize;
}

/*
 * Resets the array in the map to hold all 0 values*/
void resetKmerCountMap(KmerCountMap* kmerCountMap) {
	int* kmerCountArray = kmerCountMap->kmerCountArray;
	for (int i = 0; i < kmerCountMap->size; i++) {
		kmerCountArray[i] = 0;
	}
}

/*
 * Frees the space allocated for the kmerCountArray in the given map
 * */
void deleteKmerCountMap(KmerCountMap* kmerCountMap) {
	free(kmerCountMap->kmerCountArray);	
}

/*
 * Prints in tab format the kmers being recoreded by the kmerCountMap
 * */
void printKmers(KmerCountMap* kmerCountMap) {
	char kmer[KMER_MAX_SIZE + 1]; //space for kmer and null terminator in case size is KMER_MAX_SIZE
	int i, kmerSize = kmerCountMap->kmerSize;
	for (i = 0; i < kmerCountMap->size; i++) {
		
		indexToKmer(i, kmerSize, kmer);
		printf("%s\t", kmer);		
	
	}

}

/*
 * Prints in tab format the frequency of each kmer in the map.
 * Sums of the total kmer sequences, and divides each by that value.*/
void printFreqs(KmerCountMap* kmerCountMap) {
	double sum = 0;
	int* kmerCountArray = kmerCountMap->kmerCountArray;
	int i;
	for (i = 0; i < kmerCountMap->size; i++) {
		sum += kmerCountArray[i];
	}
	for (i = 0; i < kmerCountMap->size; i++) {
		printf("%.10lf\t", kmerCountArray[i] / sum );
	}
}

/*
 * Prints results for the kmerCountMap*/
void dumpResults(KmerCountMap* kmerCountMap) {

	int i;

	for (i = 0; i < kmerCountMap->size; i++) {
	
		printf("@index %4d: %d\t", i, kmerCountMap->kmerCountArray[i]);
	
	}

}

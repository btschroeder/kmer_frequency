#!/bin/bash

#this script compares my c version of the kmer frequency program vs the python version
#it runs each version for the same file for kmers sized 2 to size 10 by default

# Let's add some command line args

LOWER=2			#testing kmer sizes starting at 2
UPPER=10		#testing kmer sizes starting at 10
OUTPUT_DIR=outputs	#dump output of kmer tests to outputs
BENCH_DIR=benchmarks	#dump time data of kmer tests to benchmarks
PYTHON_SCRIPT=python_version/kmer_freq.py	#path to python script to run for tests (symlink)
BINARY_FILE=kmer_freq				#path to binary file testing against python script
INPUT_FILE=python_version/file_test.fasta	#path to input test file


while [ $# -gt 0 ]
do
	key=$1
	case $key in
		-l|--lower)		#specify lower kmer size bound (inclusive)
			LOWER=$2
			shift #put value in $1
			shift #now next argument is in $1
			;;
		-u|--upper)		#specify upper kmer size bound (inclusive)
			UPPER=$2
			shift
			shift
			;;
		-d|--outputs-dir)	#specify directory for outputs of tests
			OUTPUT_DIR=$2
			shift
			shift
			;;
		-b|--benchmarks-dir)	#specify directory for outputs of timing of tests
			BENCH_DIR=$2
			shift
			shift
			;;
		-p|--python-script)	#specify python script location used to run test
			PYTHON_SCRIPT=$2
			shift
			shift
			;;
		-f|--bin-file)		#specify binary file location used to run test
			BINARY_FILE=$2
			shift
			shift
			;;
		-i|--input)		#specify test file location as input for the tests
			INPUT_FILE=$2
			shift
			shift
			;;
		*)			#unspecified argument
			shift		#just ignore it
			;;
	esac
done


	

mkdir $BENCH_DIR 2> /dev/null #ignore if already exists
mkdir $OUTPUT_DIR 2> /dev/null #same

for ((size=$LOWER; size<=$UPPER; size++))
do
	/usr/bin/time -p python3 $PYTHON_SCRIPT -f $INPUT_FILE -k $size > "${OUTPUT_DIR}/python_file_test_k${size}.out" 2> "${BENCH_DIR}/python_file_test_k${size}.time"  
	/usr/bin/time -p "./${BIN_FILE}" $INPUT_FILE $size > "${OUTPUT_DIR}/c_file_test_k${size}.out" 2> "${BENCH_DIR}/c_file_test_k${size}.time"

done

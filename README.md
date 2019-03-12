
# kmer_frequency
C version of kmer_freq.py

# make
To compile the c code just run `make` in the root directory of this project.
The output is `kmer_freq`.

# run
The usage is `./kmer_freq <input_file> <kmer_size>`
Append `> OUTPUT` to the end to print to file specified by OUTPUT instead of terminal.

# benchmarks
To run benchmarks which test both the compiled C version and the Python version run `./benchmark.sh`

You can supply the following arguments to customize it:

`-l, --lower VALUE` search for kmers starting at size VALUE     (default lower bound is 2)

`-u, --upper VALUE` search for kmers up to and including VALUE  (default upper bound is 10)

`-d, --outputs-dir PATH` dump the outputs of frequency files in the directory specified by PATH

`-b, --benchmakrs-dir PATH` dumpt the outputs of the timing data to directory specified by PATH

`-p, --python-script PATH`  specify path to the python script

`-f, --bin-file PATH`       specify the path to the binary file (in this case compiled C file)

`-i, --input PATH`          specify the path to the input file used for running tests


If you do not specify any arguments this script will test kmers 2 up to 10 in size, and work correctly with this default directory structure.

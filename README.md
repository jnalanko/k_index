K-mer pseudoalignment using the r-index. Work in progress. Tested on Ubuntu 18.04.4 LTS.

Requirements: 

* CMake
* C++ compiler recent enough to have support for the <filesystem> header (e.g. g++-8).
* The sdsl-library installed at ~/lib and ~/include.

# Compiling

This project uses a modified version of Nicola Prezza's r-index implementation at https://github.com/nicolaprezza/r-index.
Compile the r-index component with:

```
# Compile r-index
cd r-index
mkdir build
cd build
cmake ..
make
cd ../..
```

Compile the rest of the code by running the following:

```
make preprocess_data 
make k_index_build 
make pseudoalign 
```

Binaries should now be at `bin`.

# Building the index

There is example data provided at the directory `genomes`. The index is built with:

```
bin/preprocess_data genomes index/concat.txt
bin/k_index_build index/concat.txt
```


This stores the index components into the directory `index`.

# Running queries against the index

There are example queries at the directory `reads`. The queries are in FASTA format. Each query is assigned
an integer identifier 0,1,2... in the same order as the reads appear in the file. To run the queries against 
the example index built above, run:

```
bin/pseudoalign index/concat.txt reads/ERR434909_1_prefix.fna
```

The output will be printed to stdout. There will be one line per query sequence. The line starts with
the identifier of the query sequence, followed by a space-separated list of integers specifying the identifiers
of the reference sequences that the query pseudoaligned with. These identifiers can be translated back to the
original filenames by looking at the file index/concat.txt.docs.names. Line i of the file is the filename of
the reference sequence with identifier i (i = 0,1,2...).

# Caveats/todos

* Currently, the value k=30 is hardcoded
* Construction is done via the suffix array and is therefore space consuming
* Documentation is poor and error handling is unhelpful

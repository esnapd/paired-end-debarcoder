# paired-end-debarcoder

paired-end-debarcoded is a python package for demultiplexing paired-end Illumina reads where custom primers have been designed to incorporate
barcodes as the first basepairs of the forward and reverse reads when using the standard Illumina sequencing primers. This package offers a few 
functions for processing the forward and reverse reads in order to demultiplex sample runs.  At its very simplest, the core function is
`demultiplexfastq` which will split the Forward/Reverse Fastq read files according to a barcode file and will output these fastqs to a 
specified directory.


## Commands In this Package

This package uses [Click](http://click.pocoo.org/6/) to make the following functions accessible on the command line.

* `fastq-concat` - Concatenate forward and reverse reads together.
* `demultiplexfastq` - Demultiplex paired-end Fastqs of PCR amplicons generated with paired-end barcodes.
* `demultiplexfasta` - Demultiplex paired-end Fastqs of PCR amplicons generated with paired-end barcodes.
* `demultiplexfasta_parallel` - Parallel version of the above.

## Recommended Installation
Intended for use with the [conda package management tool](http://conda.pydata.org/docs/index.html), and in particular with the use of their 
[environment](http://conda.pydata.org/docs/using/envs.html) files.
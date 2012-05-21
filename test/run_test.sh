#!/bin/bash

mkdir demo_run
cd demo_run

echo 'unpacking fq'
cp ../demo_1.fq.bz2 .
cp ../demo_2.fq.bz2 .
bunzip2 demo_1.fq.bz2
bunzip2 demo_2.fq.bz2

echo 'build bowtie db'
bowtie-build ../test.fasta db

echo 'running bowtie'
bowtie -I 1 -X 100000 -k 4 -S -p 4 --strata --best \
    db -1 demo_1.fq -2 demo_2.fq demo.sam


echo 'running hagfish'
hagfish_extract -S -v demo.sam 
hagfish_coverage_combine -v
hagfish_gapfinder -v -f ../test.fasta 

echo 'and plot'
hagfish_cplot2 --ymax 200 -n 2e4 contig
hagfish_blockplot -n 2e4 contig



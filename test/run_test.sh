#!/bin/bash

mkdir -p demo_run 2>/dev/null
cd demo_run

if [[ ! -f demo_1.fq ]]
then
    echo 'unpacking demo fq'
    cp ../demo_1.fq.bz2 .
    cp ../demo_2.fq.bz2 .
    bunzip2 demo_1.fq.bz2
    bunzip2 demo_2.fq.bz2
fi

if [[ ! -f alter_1.fq ]]
then
    echo 'unpacking alternative fq'
    cp ../alter_1.fq.bz2 .
    cp ../alter_2.fq.bz2 .
    bunzip2 alter_1.fq.bz2
    bunzip2 alter_2.fq.bz2
fi

if [[ ! -f db.1.ebwt ]]
then
    echo 'build bowtie db'
    bowtie-build ../test.fasta db
fi

if [[ ! -f "demo.sam" ]]
then
    echo 'running bowtie'
    bowtie -I 1 -X 100000 -k 4 -S -p 4 --strata --best \
        db -1 demo_1.fq -2 demo_2.fq demo.sam
fi

if [[ ! -f "alter.sam" ]]
then
    echo 'running bowtie'
    bowtie -I 1 -X 100000 -k 4 -S -p 4 --strata --best \
        db -1 alter_1.fq -2 alter_2.fq alter.sam
fi

echo 'running hagfish'
hagfish_extract -S -v --low 210 --high 365 demo.sam 
hagfish_extract -S -v --low 210 --high 365 alter.sam 
hagfish_coverage_combine -v
hagfish_gapfinder -v -f ../test.fasta 

echo 'and plot'
hagfish_cplot2 --ymax 400 -n 6e4 contig
hagfish_blockplot -n 6e4 contig

hagfish_compplot2 -n 6e4 -l demo -L alter



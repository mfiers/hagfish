#!/bin/bash

mkdir -p demo_run 2>/dev/null
cd demo_run

if [[ ! -f set1_1.fq ]]
then
    echo 'unpacking demo fq'
    cp ../set1_1.fq.bz2 .
    cp ../set1_2.fq.bz2 .
    bunzip2 set1_1.fq.bz2
    bunzip2 set1_2.fq.bz2
fi

if [[ ! -f set2_1.fq ]]
then
    echo 'unpacking alternative fq'
    cp ../set2_1.fq.bz2 .
    cp ../set2_2.fq.bz2 .
    bunzip2 set2_1.fq.bz2
    bunzip2 set2_2.fq.bz2
fi

if [[ ! -f db.1.ebwt ]]
then
    echo 'build bowtie db'
    bowtie-build ../test.fasta db
fi

if [[ ! -f "set1.sam" ]]
then
    echo 'running bowtie'
    bowtie -I 1 -X 100000 -k 4 -S -p 4 --strata --best \
        db -1 set1_1.fq -2 set1_2.fq set1.sam
fi

if [[ ! -f "set2.sam" ]]
then
    echo 'running bowtie'
    bowtie -I 1 -X 100000 -k 4 -S -p 4 --strata --best \
        db -1 set2_1.fq -2 set2_2.fq set2.sam
fi

echo 'running hagfish - if necessary'

[[ -d "readpairs/set1" ]] || \
    hagfish_extract -S -vv --low 210 --high 365 set1.sam 
[[ -d "readpairs/set2" ]] || \
    hagfish_extract -S -vv --low 210 --high 365 set2.sam 
[[ -d "combined" ]] || hagfish_coverage_combine -v
[[ -d "gaps" ]] || hagfish_gapfinder -v -f ../test.fasta 

echo 'and plot!'

c1='-n 6e4 --dpi 400'
c2='-n 15e3 -e 15e3 --dpi 400'
cp='--ymax 400 -S -H 400'

hagfish_cplot2 $c1 $cp -l set1  contig 
hagfish_cplot2 $c1 $cp -l set2 contig

hagfish_cplot2 $c2 $cp -l set1  contig
hagfish_cplot2 $c2 $cp -l set2 contig

hagfish_blockplot $c1 -l set1 contig
hagfish_blockplot $c1 -l set2 contig

hagfish_blockplot $c2 -l set1  contig
hagfish_blockplot $c2 -l set2  contig

hagfish_compplot2 $c1 $cp -l set1 -L set2 contig
hagfish_compplot2 $c2 $cp -l set1 -L set2 contig

hagfish_blockcompplot2 $c1 -l set1 -L set2 contig
hagfish_blockcompplot2 $c2 -l set1 -L set2 contig



#!/bin/bash

INFILE="ena_sra-run_20211126-0744.tsv"
PREFIX="ftp://ftp.sra.ebi.ac.uk/vol1/fastq"
accessions=`tail +2 $INFILE | awk '{print $1}' | sed 's/\"//g'`

chdir ENAfiles

for a in $accessions; do
	DIR1=`echo $a | cut -c1-6`
	DIR2="00`echo $a | cut -c10`"
	URL=${PREFIX}/${DIR1}/${DIR2}/$a/$a.fastq.gz
	wget $URL &
done



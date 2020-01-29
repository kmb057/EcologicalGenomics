#!/bin/bash

cd ~/Ecological_Genomics/myresults

#creating new dir to store results
mkdir fastqc

for file in /data/project_data/RS_ExomeSeq/fastq/edge_fastq/XFS*fastq.gz

do

fastqc ${file} -o fastqc/

done


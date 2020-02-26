
#!/bin/bash

cd ~/Ecological_Genomics/myresults

#creating new dir to store results
mkdir 2fastqc

for file in /data/project_data/RS_RNASeq/fastq/LOL*C*.fastq.gz

do

fastqc ${file} -o 2fastqc/

done

for file in  /data/project_data/RS_RNASeq/fastq/LOL*D*.fastq.gz

do

fastqc ${file} -o 2fastqc/

done


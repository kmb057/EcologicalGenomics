#!/bin/bash

#We'll use this wrapper to run our different scripts.

#Path to my repo

myrepo="/users/k/b/kburns/EcologicalGenomics/"

#My population
mypop="XFS"

#Directory to our cleaned and paired reads
input="/data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${mypop}"

#Directory to store the outputs on our mapping
output="/data/project_data/RS_ExomeSeq/mapping"


#Run mapping.sh
#source ./mapping.sh

#Run the post-processing steps
source ./process_bam.sh




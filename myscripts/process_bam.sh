#for f in ${output}/BWA/${mypop}*.sam

#do

#  out=${f/.sam/}
 # sambamba-0.7.1-linux-static view -S --format=bam ${f} -o ${out}.bam
  #samtools sort ${out}.bam -o ${out}.sorted.bam

#done

#Now, let's remove the PCR duplicates
#for file in ${output}/BWA/${mypop}*.sorted.bam

#do
 # f=${file/.sorted.bam/}
  #sambamba-0.7.1-linux-static markdup -r -t 1 ${file} ${f}.sorted.rmdup
#done

#Now, to finish we'll index our files
for file in ${output}/BWA/${mypop}*.sorted.rmdup.bam

do
  samtools index ${file}
done

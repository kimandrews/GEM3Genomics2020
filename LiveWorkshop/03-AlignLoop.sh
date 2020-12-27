module load bwa
module load samtools

for fastq_trim in ./02-Trim_subset/*_1.trim.sub.fastq
  do
    samplename=$(basename $fastq_trim _1.trim.sub.fastq)
    bwa mem ./Ref/ecoli_rel606.fasta \
    ./02-Trim_subset/${samplename}_1.trim.sub.fastq ./02-Trim_subset/${samplename}_2.trim.sub.fastq \
    > ./03-Align/${samplename}.aligned.sam

    samtools view -S -b ./03-Align/${samplename}.aligned.sam > ./03-Align/${samplename}.aligned.bam
    samtools sort ./03-Align/${samplename}.aligned.bam -o ./03-Align/${samplename}.aligned.sorted.bam
    samtools index ./03-Align/${samplename}.aligned.sorted.bam 
    rm ./03-Align/${samplename}.aligned.sam
    rm ./03-Align/${samplename}.aligned.bam 

done

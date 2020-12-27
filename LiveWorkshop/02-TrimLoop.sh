module load trimmomatic
for fastq in ./00-RawData/*_1.fastq.gz
do
  samplename=$(basename $fastq _1.fastq.gz)

  trimmomatic PE ./00-RawData/${samplename}_1.fastq.gz ./00-RawData/${samplename}_2.fastq.gz \
                ./02-Trim/${samplename}_1.trim.fastq.gz ./02-Trim/${samplename}_1un.trim.fastq.gz \
                ./02-Trim/${samplename}_2.trim.fastq.gz ./02-Trim/${samplename}_2un.trim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 \
                ILLUMINACLIP:/opt/modules/biology/trimmomatic/0.33/bin/adapters/NexteraPE-PE.fa:2:40:15


done

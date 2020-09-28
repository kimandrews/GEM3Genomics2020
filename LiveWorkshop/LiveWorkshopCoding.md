# GEM3 Genomics Workshop: Live Coding

We will be analyzing whole genome Illumina shotgun sequence data from the bacteria *Escherichia coli* to identify genetic variants. The data and analyses used here are modified from the publicly available Data Carpentries Genomics Workshop.

First, we need to get the raw sequence data. Go to your home directory, create a new directory for your analyses, and navigate into that directory:

```
cd ~
mkdir ./GEM3GenomicsWorkshop
cd ./GEM3GenomicsWorkshop
```
Make a directory for the raw sequence data, and navigate into the directory:
```
mkdir 00-RawData
cd 00-RawData
```

Create Symlinks (symbolic links) to the raw data, and then view the symlinks:

```
ln -s /mnt/ceph/kandrews/GEM3GenomicsWorkshopData/genomics_raw_data/* .
ls -la 
```
View one of the raw datafiles:

```
zcat SRR2584863_1.fastq.gz | less -S
```

Move back up into the GEM3GenomicsWorkshop directory:

```
cd ../
```

## Quality control
Next, we will evaluate the quality of the raw sequence data using two programs: Fastqc and Multiqc. Fastqc summarizes quality for one fastq file at a time, whereas Multiqc summarizes quality for multiple fastq files at a time.

Make a new directory to store the output of quality control analyses:

```
mkdir ./01-Quality
```

Load the Fastqc program, and run the program on all the raw sequence data files, putting the output in the 02-Qual folder:
```
module load fastqc
fastqc ./00-RawData/* -o ./01-Quality
```
Transfer one of the html output files to your local computer, and view it with a web browser.

Load a python module so you can run Multiqc, and use Multiqc to analyze and summarize all the Fastqc output files:
```
module load python/3.6.7
multiqc -i fastqc ./01-Quality -o ./01-Quality
```
Transfer one of the html output files to your local computer, and view it with a web browser.

Next, we will "clean" the sequence reads using Trimmomatic. Make a directory to store the output, then load the program:
```
mkdir ./02-Trim
module load trimmomatic
```
Run Trimmomatic for one sample:
```
trimmomatic PE ./00-RawData/SRR2589044_1.fastq.gz ./00-RawData/SRR2589044_2.fastq.gz \
                ./02-Trim/SRR2589044_1.trim.fastq.gz ./02-Trim/SRR2589044_1un.trim.fastq.gz \
                ./02-Trim/SRR2589044_2.trim.fastq.gz ./02-Trim/SRR2589044_2un.trim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 \
                ILLUMINACLIP:/opt/modules/biology/trimmomatic/0.33/bin/adapters/NexteraPE-PE.fa:2:40:15 
```
## Align to a reference genome
Next, we will align the cleaned sequence reads to a reference genome for one sample using the program BWA.

First, create a directory and a symlink for the reference genome:
```
mkdir ./Ref
ln -s /mnt/ceph/kandrews/GEM3GenomicsWorkshopData/genomics_Ref/ecoli_rel606.fasta ./Ref
ls -la ./Ref
```
Load BWA, and index the reference genome:
```
module load bwa
bwa index ./Ref/ecoli_rel606.fasta
```
Next we will align the cleaned sequence reads to the reference genome. We will be using a subsampled dataset to make the analysis go faster. The subsampled dataset consists of a smaller number of cleaned reads per sample.

Make a directory for the subsampled dataset, create symlinks to the dataset, and view the symlinks:
```
mkdir 02-Trim_subset
ln -s /mnt/ceph/kandrews/GEM3GenomicsWorkshopData/genomics_trim_subset/*  ./02-Trim_subset
ls -la ./02-Trim_subset
```
Make a directory to store the results of genome alignment, and align the reads from one sample to the reference genome:
```
mkdir ./03-Align
bwa mem ./Ref/ecoli_rel606.fasta ./02-Trim_subset/SRR2584866_1.trim.sub.fastq ./02-Trim_subset/SRR2584866_2.trim.sub.fastq > ./03-Align/SRR2584866.aligned.sam
```
View the output, which is a sam file:
```
less -S ./03-Align/SRR2584866.aligned.sam
```
Convert the sam file to a bam file, and then sort and index the bam file, using the program samtools:
```
module load samtools

samtools view -S -b ./03-Align/SRR2584866.aligned.sam > ./03-Align/SRR2584866.aligned.bam

samtools sort ./03-Align/SRR2584866.aligned.bam -o ./03-Align/SRR2584866.aligned.sorted.bam

samtools index ./03-Align/SRR2584866.aligned.sorted.bam 
```
Now, remove the sam file and the unsorted bam file, to save hard drive space:
```
rm ./03-Align/SRR2584866.aligned.sam
rm ./03-Align/SRR2584866.aligned.bam 
```

Learn more about the bam file:
```
samtools flagstat ./03-Align/SRR2584866.aligned.sorted.bam 
```
## Variant Calling (Genotyping)

Identify and genotype variants for one sample using bcftools:

```
module load bcftools
mkdir ./04-Genotype

bcftools mpileup -O b -f ./Ref/ecoli_rel606.fasta \
./03-Align/SRR2584866.aligned.sorted.bam -o ./04-Genotype/SRR2584866_raw.bcf 

bcftools call --ploidy 1 -m -v ./04-Genotype/SRR2584866_raw.bcf  -o ./04-Genotype/SRR2584866_variants.vcf 
```
View the output vcf file:
```
less -S ./04-Genotype/SRR2584866_variants.vcf
```

## Analyze all the samples using shell scripts

Use shell scripts to analyze all of the samples efficiently.

Copy the scripts to your account:
```
cp /mnt/ceph/kandrews/GEM3GenomicsWorkshopScripts/02-TrimLoop.sh ~/GEM3GenomicsWorkshop

cp /mnt/ceph/kandrews/GEM3GenomicsWorkshopScripts/03-AlignLoop.sh ~/GEM3GenomicsWorkshop
```

Run Trimmomatic for all samples:
```
bash 02-TrimLoop.sh
```
Align all samples to the reference genome, using the subsetted, trimmed datafiles so the analysis will go faster:
```
bash 03-AlignLoop.sh
```
Identify and gentype variants for all subsetted samples:
```
bcftools mpileup -O b -f ./Ref/ecoli_rel606.fasta \
./03-Align/*.aligned.sorted.bam -o ./04-Genotype/all_raw.bcf -a DP

bcftools call --ploidy 1 -m -v ./04-Genotype/all_raw.bcf  -o ./04-Genotype/all_variants.vcf 
```
Filter the variants using VCFTools.

Remove indels:
```
module load vcftools

vcftools --vcf ./04-Genotype/all_variants.vcf --out ./04-Genotype/all_variants_SNPs --remove-indels --recode
```
Remove any loci with mean depth of coverage (across samples) below 5:
```
vcftools --vcf ./04-Genotype/all_variants_SNPs.recode.vcf --out ./04-Genotype/all_variants_SNPs_dp5 --min-meanDP 5 --recode
```








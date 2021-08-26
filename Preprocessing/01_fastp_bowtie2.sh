#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
set -e
set -u
set -o pipefail
set -o nounset
set -o errexit

#Loading all needed modules
##This script was adapted from Roxana Shen's MEDIPS_PE_pipe.sh script by SLW
module load parallel
#module load fastp/0.19.4
module load fastqc/0.11.5
module load bowtie2/2.3.5
module load samtools/1.3.1
module load picard/2.10.9  
module load bedtools/2.27.1 
module load igenome-human/hg38

#mkdir data/2019_TrimmedControls/ #To store the trimmed files in later
##Make sure this is storing data sile in the working directory and not on mordor

#Checking to see what the quality of the unaligned reads are. 
parallel -j 6 "fastqc -o data/2019_SpikeInControl_data --noextract -f fastq {1}" ::: data/2019_SpikeInControl_data/*.gz

##fastp- trimming the adapters and appending UMI to sample same
## -i input of read 1
## -I input of read 2
## -o output of read 1
## -O output of read 2
## --linksm read 1 and read 2
## fastp automatically gzips outputs
parallel -j 6 --link "fastp -i {1} -I {2} -o data/2019_TrimmedControls/trimmed.{1/} -O data/2019_TrimmedControls/trimmed.{2/} --umi --umi_loc=per_read --umi_len=5 --adapter_sequence=AATGATACGGCGACCACCGAGATCTACACATATGCGCACACTCTTTCCCTACACGAC --adapter_sequence_r2=CAAGCAGAAGACGGCATACGAGATACGATCAGGTGACTGGAGTTCAGACGTGT -j data/2019_TrimmedControls/{1/.}.{2/.}.json -h data/2019_TrimmedControls/{1/.}.{2/.}.html" ::: data/2019_SpikeInControl_data/*R1*fastq.gz ::: data/2019_SpikeInControl_data/*R2*fastq.gz

##Align reads to reference genome
### Create my own index to align to the synthetic fragment sequences located in the path where bowtie2 is located
#bowtie2-build -f SyntheticDNA_seq.fasta  SyntheticDNA_idx

### -f specifies the reference genome to be used, fastq file
### -p specifies the number of threads to be used in parallelization
### --minins specifies the minimum valid fragment length required for paired end sequence- I have removed a min, as I have 80bp control fragments I want to pick up.
### --maxins specifies the maximum valid fragment length required for paired end sequence
### -q specifies the input fastq read files to be aligned, -1 specified first mate pair, -2 specifies second mate pair
### -S writes the alignment results to a sam file as output
### --no-unal suppresses SAM files for reads that failed to align
### --un will write unaligned reads to a fastq file

#mkdir data/2019_ControlAlignment/ #To store both aligned and unaligned files

##Aligns sequences to control fragments
parallel -j 6 --link "bowtie2 --local -x ./SyntheticDNA_idx -p 6 --minins 80 --maxins 320 --no-unal --un-conc-gz data/2019_ControlAlignment/unalignedtospike.{1/} -q -1 {1} -2 {2} -S data/2019_ControlAlignment/aligned.{1/.}{2/.}.sam 2> data/2019_ControlAlignment/{1/.}_bowtie2.log"  ::: data/2019_TrimmedControls/*R1*.fastq.gz ::: data/2019_TrimmedControls/*R2*.fastq.gz  

##Now align to remaining sequencing to the human genome 
parallel -j 6 --link "bowtie2 --local -x $BOWTIE2INDEX -p 6 --minins 80 --maxins 320 -q -1 {1} -2 {2} -S data/2019_ControlAlignment/human.aligned.{1/.}{2/.}.sam 2> data/2019_ControlAlignment/{1/.}_bowtie2.log" ::: data/2019_ControlAlignment/unalignedtospike*fastq.1.gz ::: data/2019_ControlAlignment/unalignedtospike*fastq.2.gz

#Convert .sam to .bam files
parallel -j 6 "samtools view -bSh {1} > data/2019_ControlAlignment/{1/.}.bam" ::: data/2019_ControlAlignment/*.sam 
## I have not sorted or indexed these files as they are synthetic DNA fragments not aligning to human genome, therefore have no coordinates. Because this is a small file, indexing is not needed.

##Mapping info
parallel -j 6 "samtools flagstat {1} 2> mappinginfo_{1/.}.log" ::: data/2019_ControlAlignment/*.bam

##Filtering reads to only those that were mater properly and counting the number of overlap between reads tha reference sequence (this should be read counts)
parallel -j 6 "samtools view -bf 0x2 {1} > data/2019_ControlAlignment/filtered_{1/.}.bam" ::: data/2019_ControlAlignment/*.bam 

#sort bam
### This something has a parallelization issue, so watch that this is correct
parallel -j 12 "samtools sort {1} {1.}_sorted" ::: /mnt/work1/users/home2/wilsons/Projects/2018_PTB/data/2019_ControlAlignment/*.bam 

#index bam files
parallel -j 12 "samtools index {1}" ::: /mnt/work1/users/home2/wilsons/Projects/2018_PTB/data/2019_ControlAlignment/*sorted.bam 

##UMI dedup
parallel -j 12 "umi_tools dedup -I {1} --paired -S {1.}_dedup.bam --output-stats={1.}_stats.tsv --mapping-quality=20 --unpaired-reads=discard --chimeric-pairs=discard --log={1.}.log" ::: /cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/*sorted.bam


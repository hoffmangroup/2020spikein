#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
set -e
set -u
set -o pipefail
set -o nounset
set -o errexit

#Loading all needed modules
module load parallel
module load samtools
module load bedtools

#convert deduplicated bam to bedpe
parallel -j 12 "samtools sort -n {1} -o {1.}_sort.bam" ::: /cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/*dedup.bam

#convert deduplicated bam to bedpe
parallel -j 12 "bedtools bamtobed -bedpe -i {1}" ::: /cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/*dedup.bam

#cut hg38 aligned bed file to only to read 1 start and read 2 end
parallel -j 6 "cut -f1-2,6,7,11- {1} > {1.}_cut.bed" ::: /cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/filtered_human*.bed

#calculate GC content, fragment length, and CpG number
## As some CG are capital and some not, I do this twice

parallel -j 6 "bedtools nuc -fi /cluster/projects/hoffmangroup/data_samanthawilson/2020_PTB_secondtri/hg38.fa -bed {1} -pattern CG > {1.}_CG.bed " ::: /cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/*_cut.bed

parallel -j 6 "bedtools nuc -fi /cluster/projects/hoffmangroup/data_samanthawilson/2020_PTB_secondtri/hg38.fa -bed {1} -pattern CG > {1.}_cg.bed " ::: /cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/*_cut.bed

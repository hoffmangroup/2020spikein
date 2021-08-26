#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
set -e
set -u
set -o pipefail
set -o nounset
set -o errexit

##Find the overlaps between our fragments and each of the 300bp windows in genom
### Retrun a file with the fragment listed, the concentration value for each sample, the window it overlaps and the number of basepairs that it overlaps

module load bedtools
module load parallel


parallel -j 6 "bedtools intersect -wao -a {1} -b /cluster/home/wilsons/Projects/2020_Control_Project/2020_BatchAnalysis/2020_hg38_300bpwindows.bed > {1.}_hg38_intersect.bed" ::: /cluster/projects/hoffmangroup/data_samanthawilson/*pmol.bed

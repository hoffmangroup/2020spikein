#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
set -e
set -u
set -o pipefail
set -o nounset
set -o errexit

module load parallel

##Check that the first 3 rows match between the files
parallel -j 6 "awk '{k=$1 FS $2 FS $3 FS $4} NR==FNR{a[k]; next} !(k in a)' {1} {2}" ::: ~/Projects/2018_PTB/data/2019_ControlAlignment/*CG.bed ::: ~/Projects/2018_PTB/data/2019_ControlAlignment/*cg.bed
##event not found, they should match

##Paste last column of cg file to the CG file
parallel -j 1 "cut -f14-14 {1} > {1.}_tmp.bed" ::: ~/Projects/2018_PTB/data/2019_ControlAlignment/*cg.bed

parallel -j 1 "paste {1} {2} > {1.}_all.bed" ::: ~/Projects/2018_PTB/data/2019_ControlAlignment/*CG.bed ::: ~/Projects/2018_PTB/data/2019_ControlAlignment/*tmp.bed

##Make sure to delete tmp files after



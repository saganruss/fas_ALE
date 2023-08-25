#!/bin/bash

# load the snippy env
source ~/miniconda3/etc/profile.d/conda.sh
# mamba init
conda activate snippy_env

#set path variables
genome='/proj/omics4tb2/sruss/fas_ALE/data/genome/genomic.gbff'
out='/proj/omics4tb2/sruss/fas_ALE/data/snippy_results/'

# run snippy on each sample seperately
snippy --outdir ${out}${3}_snps --ref ${genome} --R1 ${1} --R2 ${2}

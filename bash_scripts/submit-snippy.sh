#!/bin/bash
#SBATCH --array=1-20
#SBATCH -J snippy
#SBATCH -o "%A"."%a".out
#SBATCH -e "%A"."%a".out
#SBATCH --mem-per-cpu=1200
#SBATCH -c 8

# Job array for variant calling -- snippy aligns the read pairs to the reference genome and then calls SNPs based on this alignment

# load the snippy env
source ~/miniconda3/etc/profile.d/conda.sh
# mamba init
conda activate snippy_env

# set variable for parsing SLURM array and sample list
samplesheet="/proj/omics4tb2/sruss/fas_ALE/data/directories.tsv"

r1=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $1}'`
r2=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $2}'`
name=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $3}'`

#run the snippy script for each sample separately
echo "processing sample" ${name} "with reads" ${r1} "and" ${r2}
srun /proj/omics4tb2/sruss/fas_ALE/code/bash_scripts/run-snippy.sh ${r1} ${r2} ${name}

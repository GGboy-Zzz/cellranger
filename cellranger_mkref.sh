#!/bin/bash
#SBATCH --job-name=mkref
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # specify number of CPUs to use here
#SBATCH --mem-per-cpu=10gb

module load cellranger/7.1.0
cellranger mkref  --genome=mm39_hsv1 \
--fasta=/home/zhuyong/cellranger/00ref/GCF_000001635.27_GRCm39_genomic_HSV1.fna \
--genes=/home/zhuyong/cellranger/00ref/mm39_hsv1_filtered.gtf


#!/bin/bash
#SBATCH --job-name=OHSV3
#SBATCH --partition cpu
#SBATCH -o /home/zhuyong/cellranger/code/logs/%j_%x.log 
#SBATCH -e /home/zhuyong/cellranger/code/logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30 # specify number of CPUs to use here
#SBATCH --mem-per-cpu=10gb
module load cellranger/7.1.0

cellranger count --id=OHSV3  \
                 --transcriptome=/home/zhuyong/cellranger/00ref/mm39_hsv1  \
                 --fastqs=/data/share_for_user/1207_zhangjunjLab/DZOE2023051540-b1/OHSV3

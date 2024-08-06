###### Cellranger tutorial 
# ZY
# 27 Oct 2023
# refer to: https://www.10xgenomics.com/support/single-cell-gene-expression

##### Linux script -------

### Build a Custom Reference (cellranger mkref)
# find the input files, filter the GTF, remove entries of non-ployA transcripts
cellranger mkgtf /home/zhuyong/cellranger/00ref/mm39.gtf  /home/zhuyong/cellranger/00ref/mm39_filtered.gtf   \
                   --attribute=gene_biotype:protein_coding \
                   --attribute=gene_biotype:lncRNA \
                   --attribute=gene_biotype:antisense \
                   --attribute=gene_biotype:IG_LV_gene \
                   --attribute=gene_biotype:IG_V_gene \
                   --attribute=gene_biotype:IG_V_pseudogene \
                   --attribute=gene_biotype:IG_D_gene \
                   --attribute=gene_biotype:IG_J_gene \
                   --attribute=gene_biotype:IG_J_pseudogene \
                   --attribute=gene_biotype:IG_C_gene \
                   --attribute=gene_biotype:IG_C_pseudogene \
                   --attribute=gene_biotype:TR_V_gene \
                   --attribute=gene_biotype:TR_V_pseudogene \
                   --attribute=gene_biotype:TR_D_gene \
                   --attribute=gene_biotype:TR_J_gene \
                   --attribute=gene_biotype:TR_J_pseudogene \
                   --attribute=gene_biotype:TR_C_gene

# Add exogenous sequences to a custon reference
cp mm39_filtered.gtf mm39_hsv1_filtered.gtf
cp mm39.fa mm39_hsv1.fa

#find the number of bases
cat hsv1.fa | grep -v "^>" | tr -d "\n" | wc -c
#make a custom GTF
echo -e 'GLU734771.1\tunknown\texon\t1\t152151\t.\t+\t.\tgene_id "hsv1"; transcript_id "hsv1"; gene_name "hsv1"; gene_biotype "protein_coding";' > hsv1.gtf
cat hsv1.gtf

cat hsv1.gtf >> mm39_hsv1_filtered.gtf
cat hsv1.fa >> mm39_hsv1.fa

tail mm39_hsv1_filtered.gtf
tail mm39_hsv1.fa
#To confirm that the GFP entry was added to the FASTA file
grep ">" mm39_hsv1.fa

cellranger mkref \
--genome=mm39_hsv1 \
--fasta=mm39_hsv1.fa \
--genes=mm39_hsv1_filtered.gtf


### Run cellranger count
cellranger count --id= \
                 --transcriptom= ref \
                 --fastqs= fastq \
                 --samples=  \      #要选择的fastq文件的前缀，optional
                 --localcores=4 \   #optional
                 --localmem=64      #optional

### Slurm script
#!/bin/bash
#SBATCH --job-name=
#SBATCH -o logs/%j_%x.log
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # specify number of CPUs to use here
#SBATCH --mem-per-cpu=10gb
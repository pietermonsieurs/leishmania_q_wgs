#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=4
#SBATCH --time=70:30:00
#SBATCH -A ap_itg_mpu

# load necessary software settings for GATK 
module load Java
module load GATK
module load BioTools

## hard coded example
# bam_file=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_plosntds/VL9.mapq30.removedups.proper_paired.bam
ref_genome=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/refgenomes/Leishmania_donovani_16Nov2015beta.fa

## parameter settings
threads=4

## index the bam file - if needed
samtools index ${bam_file}

## GATK using the GVCF approach
vcf_file=$(echo $bam_file | sed -e 's/.mapq30.removedups.proper_paired.bam$/.gatk.g.vcf/')
gatk HaplotypeCaller -R $ref_genome -I $bam_file -O $vcf_file -ERC GVCF --native-pair-hmm-threads $threads


## run for all bam files in the Somalileish folder
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/results/bwa/
# for bam_file in ${PWD}/*paired.bam; do sbatch --export=bam_file=${bam_file} /user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/bin/GATK_refgenome.sh; done
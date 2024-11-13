#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=2
#SBATCH --time=0:40:00
#SBATCH -A ap_itg_mpu

## load modules
module load BWA
module load BioTools
module load Java

## run deeptools to look for CNVs. This should be done on the 
## bam files which are not corrected for duplicate reads etc as
## the mapping quality will often be 0 due to multimapping reads
cov_file=${bam_file/.bam/.cov}
samtools depth -a ${bam_file} > ${cov_file}
# gzip ${cov_file}

## run for all .bam files in the folder of Allison 
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/results/
# for sample in {001..036}; do sbatch --export=bam_file=${PWD}/${sample}.bam /user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/bin/cnv_samtools_depth.slurm; done
# for sample in {001..012}; do sbatch --export=bam_file=${PWD}/${sample}_Yeti.bam /user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/bin/cnv_samtools_depth.slurm; done




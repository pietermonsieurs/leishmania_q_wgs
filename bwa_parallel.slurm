#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=10:00:00
#SBATCH -A ap_itg_mpu

## load modules
module load BWA
module load BioTools
module load Java

## define variables
export THREADS=28
export SEED=50

## BWA_DIR
BWA_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/results/bwa/

## REF GENOME
export REF_GENOME=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/refgenomes/Leishmania_donovani_16Nov2015beta.fa

## SOFTWARE PARAMS
export PICARD_JAR=/data/antwerpen/205/vsc20587/software/picard/picard.jar
export PICARD_BIN="java -jar ${PICARD_JAR}"
export mapq_cutoff=30


## create the necessary filenames based on te input information
# fastq_file_R1=/user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/data/fastq/028_R1.fastq.gz
echo "before substituions: ${fastq_file_R1}"
fastq_file_R2=${fastq_file_R1/_R1.fastq.gz/_R2.fastq.gz}

## file prefix - depening on naming convention
file_prefix=${fastq_file_R1%_R1.fastq.gz}
file_prefix=${file_prefix##*/}
bam_file_prefix=${BWA_DIR}/${file_prefix}

echo $fastq_file_R1
echo $fastq_file_R2
echo $file_prefix
echo $bam_file_prefix

### run BWA with increased seed 
bwa mem \
    -R "@RG\tID:${file_prefix}\tSM:${file_prefix}\tPL:ILLUMINA" \
    -k $SEED \
    -t $THREADS \
    $REF_GENOME \
    $fastq_file_R1 $fastq_file_R2 | \
    samtools sort -@$THREADS -o ${bam_file_prefix}.bam
samtools index ${bam_file_prefix}.bam
samtools flagstat -@$THREADS ${bam_file_prefix}.bam > ${bam_file_prefix}.flagstat


### selection on mapping quality ###
samtools view -@ $THREADS -bq $mapq_cutoff ${bam_file_prefix}.bam > ${bam_file_prefix}.mapq${mapq_cutoff}.bam
samtools index ${bam_file_prefix}.mapq${mapq_cutoff}.bam
samtools flagstat ${bam_file_prefix}.mapq${mapq_cutoff}.bam > ${bam_file_prefix}.mapq${mapq_cutoff}.flagstat


### remove duplicate reads ### 
$PICARD_BIN MarkDuplicates REMOVE_DUPLICATES=true I=${bam_file_prefix}.mapq${mapq_cutoff}.bam O=${bam_file_prefix}.mapq${mapq_cutoff}.removedups.bam M=${bam_file_prefix}.mapq${mapq_cutoff}.markdups_metrics.txt
samtools index ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.bam
samtools flagstat ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.bam > ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.flagstat


### select only proper paired reads ###
samtools view -@$THREADS -bf 0x2 ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.bam > ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.proper_paired.bam 
samtools index ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.proper_paired.bam 
samtools flagstat ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.proper_paired.bam > ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.proper_paired.flagstat




#### run using different sbatch commands ####

## run for the first batch of samples of the WGS data of Allsion
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/data/fastq/
# for fastq_file_R1 in  *_R1.fastq.gz; do sbatch --export=fastq_file_R1=${PWD}/${fastq_file_R1} /user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/bin/bwa_parallel.slurm; done

## run for the second batch of samples - Yeti samples - of Allison
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/data/fastq/
# for fastq_file_R1 in *_Yeti_R1.fastq.gz; do sbatch --export=fastq_file_R1=${PWD}/${fastq_file_R1} /user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/bin/bwa_parallel.slurm; done

#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=10:00:00
#SBATCH -A ap_itg_tryps

## load modules
module load Python/3.12.3-GCCcore-13.3.0
module load BWA
module load BioTools
module load BioTools-Python
module load Java

## define variables
export threads=28
export seed=50

## BWA_DIR
bwa_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_10X_Q/results/bulk/bwa/

## REF GENOME
ref_genome=/user/antwerpen/205/vsc20587/scratch/leishmania_10X_Q/data/refgenome/Leishmania_donovani_16Nov2015beta.fa
gtf_file=/user/antwerpen/205/vsc20587/scratch/leishmania_10X_Q/data/refgenome/Leishmania_donovani_versionGenomeScan.gtf

### index reference genomes
## indexing already performed in previous analysis
# bwa index $REF_GENOME

echo "before substituions: ${fastq_file_R1}"

## for in-house data
fastq_file_R2=${fastq_file_R1/_R1.fastq.gz/_R2.fastq.gz}
file_prefix=${fastq_file_R1%_R1.fastq.gz}
file_prefix=${file_prefix##*/}
bam_file_prefix=${bwa_dir}/${file_prefix}

echo $fastq_file_R1
echo $fastq_file_R2
echo $file_prefix
echo $bam_file_prefix

### run BWA with increased seed 
# bwa mem \
#     -R "@RG\tID:${file_prefix}\tSM:${file_prefix}\tPL:ILLUMINA" \
#     -k $seed \
#     -t $threads \
#     $ref_genome \
#     $fastq_file_R1 $fastq_file_R2 | \
#     samtools sort -@$threads -o ${bam_file_prefix}.bam
# samtools index ${bam_file_prefix}.bam
# samtools flagstat -@$threads ${bam_file_prefix}.bam > ${bam_file_prefix}.flagstat

## run htseq-count for all samples. Make sure to put the -t flag as the 
## gtf file does not contain exon features, only CDS features  
htseq-count -t exon -f bam -r pos -t CDS --stranded=no ${bam_file_prefix}.bam ${gtf_file} > ${bam_file_prefix}.htseq.strand_no.csv



#### run using different sbatch commands ####
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_10X_Q/data/bulk/
# for fastq_file_R1 in ${PWD}/*_R1.fastq.gz; do sbatch --export=fastq_file_R1=${fastq_file_R1} /user/antwerpen/205/vsc20587/scratch/leishmania_10X_Q/bin/bulk_bwa_htseq.slurm; done
#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=10:00:00

## load modules
module load BWA
module load BioTools
module load Java

## define variables
export THREADS=28
export SEED=50

## BWA_DIR
BWA_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/lbraz_susl_swga/cnv/


## REF GENOME
export REF_GENOME=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/refgenomes/TriTrypDB-57_LbraziliensisMHOMBR75M2904_2019_Genome.fasta

## create the necessary filenames based on te input information
# fastq_file_R1=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/20210622/HCCHFDSX2_104509-001-003_GCCAAGAC_L001_R1.fastq.gz
# fastq_file_R1=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/20210622/HCCHFDSX2_104509-001-001_AACGTGAT_L002_R1.fastq.gz
# fastq_file_R1=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/elife_fastq/ERR895175_1.fastq.gz

echo "before substituions: ${fastq_file_R1}"
## depending on which sequencing file, you have to choose either R2 or R3
# fastq_file_R2=${fastq_file_R1/_R1.fastq.gz/_R2.fastq.gz}
# fastq_file_R2=${fastq_file_R1/_R1.fastq.gz/_R3.fastq.gz}
fastq_file_R2=${fastq_file_R1/_1.fastq.gz/_2.fastq.gz}
# fastq_file_R2=${fastq_file_R1/_R1_001.fastq.gz/_R2_001.fastq.gz}

## file prefix - depening on naming convention
# file_prefix=${fastq_file_R1%_L00*_R1.fastq.gz}
# file_prefix=${fastq_file_R1%_R1.fastq.gz}
file_prefix=${fastq_file_R1%_1.fastq.gz}
# file_prefix=${fastq_file_R1%_R1_001.fastq.gz}

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


## run deeptools to look for CNVs. This should be done on the 
## bam files which are not corrected for duplicate reads etc as
## the mapping quality will often be 0 due to multimapping reads
samtools depth -a ${bam_file_prefix}.bam > ${bam_file_prefix}.cov
gzip ${bam_file_prefix}.cov


## run using different sbatch commands. Run for all fastq files
## that are produced for sWGA and SuSL
## fastq R1 files from the excel file Lbraziliensis.direct_sequencing.xlsx
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/lbraz_susl_swga/cnv
# while read -r fastq_file_R1; do sbatch --export=fastq_file_R1=$fastq_file_R1 /user/antwerpen/205/vsc20587/scratch/leishmania_susl/bin/lbraz_CNV_bwa_and_samtools.slurm;done < fastq_files_R1.csv

## run for all WGS data of Senne's publication
# cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/VIANNIA_GENOME_PROJECTS/FASTQ
# ls $PWD/*_1.fastq.gz > /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/lbraz_susl_swga/cnv/fastq_files_wgs_R1.csv
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/lbraz_susl_swga/cnv/
# while read -r fastq_file_R1; do sbatch --export=fastq_file_R1=$fastq_file_R1 /user/antwerpen/205/vsc20587/scratch/leishmania_susl/bin/lbraz_CNV_bwa_and_samtools.slurm;done < fastq_files_wgs_R1.csv





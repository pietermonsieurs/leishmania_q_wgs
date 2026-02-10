#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=14
#SBATCH --time=20:00:00
#SBATCH -A ap_itg_mpu


module load BioTools
module load Java

# parameter settings
out_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_10X_Q/results/bulk/fastQC/

## fastQC settings
fastqc_bin=/user/antwerpen/205/vsc20587/data/software/FastQC/fastqc
threads=14

# run FastQC based on the worker moduel
$fastqc_bin --outdir ${out_dir} --threads $threads $fastq_file


## run on the input directory
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_10X_Q/data/bulk/
# for fastq_file in *.fastq.gz; do sbatch --export=fastq_file=${fastq_file} /user/antwerpen/205/vsc20587/scratch/leishmania_10X_Q/bin/bulk_input_fastQC.slurm; done




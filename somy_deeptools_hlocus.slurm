#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=14
#SBATCH --time=04:00:00
#SBATCH -A ap_itg_mpu


## installation of deeptools: apparently the -I option was needed to force re-installation of some 
## libraries. Not clear whether this was a deeptools issue, or an issue with a previous 
## installation process
# module load Python/3
# export PYTHONPATH="/user/antwerpen/205/vsc20587/data/software/python_lib/lib/python3.8/site-packages/:${PYTHONPATH}"
# pip3 install -I  --prefix=/user/antwerpen/205/vsc20587/data/software/python_lib/ deeptools


## specify the deeptools binary and export the PythonPath so it can be run
## with a local binfile
module load Python/3.10.4-GCCcore-11.3.0
module load atools
module load BioTools
export PYTHONPATH="/user/antwerpen/205/vsc20587/data/software/python_lib/lib/python3.10/site-packages/:${PYTHONPATH}"
bamcoverage_bin=/user/antwerpen/205/vsc20587/data/software/python_lib/bin/bamCoverage

## load the bam-files using the atools option
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/results/bwa/bam_files.csv)

## normalization procedure: CPM should work, as this normalizes for the difference
## in sequencing depth. It does not take into account the feature length (in case of RNA-seq)
## but as the main aim is to check the sequencing depth per bin, and the binsize is a constant
## factor, CPM should be sufficient, and is easy interpretable.
## ==> --normalizeUsing CPM

## test file for bam
# bam_file=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_20210622/HCCHFDSX2_104509-001-001_AACGTGAT.mapq30.removedups.proper_paired.bam
echo bam_file
echo $bam_file
samtools index $bam_file

## parameters for running & creating of name for output
threads=14
bin_size=100
output_file=${bam_file/.bam/.deeptools.Ld23.bin${bin_size}.bedgraph}
# output_file="${output_file##*/}"

# output_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_elife/hlocus/
# output_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_aethiopica_skinslit/
# output_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_20230116/
# output_file=${output_dir}/${output_file}

$bamcoverage_bin \
    --bam $bam_file \
    --outFileName $output_file \
    --outFileFormat bedgraph \
    --binSize $bin_size \
    --normalizeUsing CPM \
    --numberOfProcessors $threads \
    --region Ld23
    # --skipNonCoveredRegions 

## run for all bam-files
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/results/bwa/
# ls $PWD/*proper_paired.bam  > bam_files.csv 
# vi bam_files.csv ## add header "bam_file"
# echo $PWD/bam_files.csv
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/results/bwa/bam_files.csv) /user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/bin/somy_deeptools_hlocus.slurm

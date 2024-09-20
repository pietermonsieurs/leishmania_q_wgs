#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=14
#SBATCH --time=04:00:00


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
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_aethiopica_skinslit/bam_files.csv)

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
output_file=${bam_file/.bam/.deeptools.bin${bin_size}.bedgraph}
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
    --numberOfProcessors $threads
    # --skipNonCoveredRegions 

## run for all bam-files
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_2022
# ls $PWD/*proper_paired.bam  > bam_files.csv 
# vi bam_files.csv ## add header "bam_file"
# echo $PWD/bam_files.csv
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_2022/bam_files.csv) /user/antwerpen/205/vsc20587/scratch/leishmania_susl/bin/check_coverage_deeptools_bamCoverage.slurm

## run for erroneous bam_file clinical Ltropica with SuSL
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_2022
# vi bam_files_crashed.csv ## add the two big files FJ20_09 and FJ20_13
# echo /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_2022/bam_files_crashed.csv
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_2022/bam_files_crashed.csv) /user/antwerpen/205/vsc20587/scratch/leishmania_susl/bin/check_coverage_deeptools_bamCoverage.slurm


## run for all bam-file of tropica (WGS data!)
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_tropica_clinical/
# ls $PWD/*proper_paired.bam  > bam_files.csv 
# vi bam_files.csv ## add header "bam_file"
# echo $PWD/bam_files.csv
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_tropica_clinical/bam_files.csv) /user/antwerpen/205/vsc20587/scratch/leishmania_susl/bin/check_coverage_deeptools_bamCoverage.slurm

## run for the previous sammples with removal of methylated DNA via digestion
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_2022_digestion/
# ls $PWD/*proper_paired.bam  > bam_files.csv 
# vi bam_files.csv ## add header "bam_file"
# echo $PWD/bam_files.csv
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_2022_digestion/bam_files.csv) /user/antwerpen/205/vsc20587/scratch/leishmania_susl/bin/check_coverage_deeptools_bamCoverage.slurm


## run for the samples where different primer sets are combined
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_2022_combined/
# ls $PWD/*proper_paired.bam  > bam_files.csv 
# vi bam_files.csv ## add header "bam_file"
# echo $PWD/bam_files.csv
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_2022_combined/bam_files.csv) /user/antwerpen/205/vsc20587/scratch/leishmania_susl/bin/check_coverage_deeptools_bamCoverage.slurm

## run for the tropica samples - combining the WGS and SureSelect data
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/tropica
# ls $PWD/*proper_paired.bam  > bam_files.csv 
# vi bam_files.csv ## add header "bam_file"
# echo $PWD/bam_files.csv
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/tropica/bam_files.csv) /user/antwerpen/205/vsc20587/scratch/leishmania_susl/bin/check_coverage_deeptools_bamCoverage.slurm

## run for the clinical isolates of february 2023
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_20230116/
# ls $PWD/*proper_paired.bam  > bam_files.csv 
# vi bam_files.csv ## add header "bam_file"
# echo $PWD/bam_files.csv
# module load atools
# sbatch --array $(arange --no_sniffer --data //user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_20230116/bam_files.csv) /user/antwerpen/205/vsc20587/scratch/leishmania_susl/bin/check_coverage_deeptools_bamCoverage.slurm


## run for the analysis of the H-locus
## create input file using the code in outbreak_Hlocus_prepare_bamfiles.sh
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_elife/hlocus/bam_files.csv) /user/antwerpen/205/vsc20587/scratch/leishmania_susl/bin/check_coverage_deeptools_bamCoverage.slurm


## run for the analysis of the Laethiopica skinslit and biopsie samples
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_aethiopica_skinslit/
# ls $PWD/*proper_paired.bam  > bam_files.csv 
# vi bam_files.csv ## add header "bam_file"
# echo $PWD/bam_files.csv
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_aethiopica_skinslit/bam_files.csv) /user/antwerpen/205/vsc20587/scratch/leishmania_susl/bin/check_coverage_deeptools_bamCoverage.slurm


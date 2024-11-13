#!/usr/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --time=10:00:00

## load the required modules
module load SciPy-bundle/2023.07-gfbf-2023a ## numpy
module load IPython/8.5.0-GCCcore-11.3.0 ## tornado
module load matplotlib ## matplotlib

## create an extended bed file. It also needs to contain the strand + 
## the gene ID
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/refgenomes/
# awk '!/^#/ && ($3 == "CDS" || $3 == "rRNA") {split($9, id, "="); split(id[2], gene_id, ";"); print $1 "\t" $4 "\t" $5 "\t" $7 "\t" gene_id[1]}' TriTrypDB-57_LbraziliensisMHOMBR75M2904_2019.gff > TriTrypDB-57_LbraziliensisMHOMBR75M2904_2019.genes.ext.bed

## do zipping in case still needed and rename the
## .cov file with an extension .gz
gzip ${cov_file}
cov_file=${cov_file}.gz

## load some test files
# cov_file=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/lbraz_susl_swga/cnv/104989-001-009.cov.gz
bed_file=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/refgenomes/TriTrypDB-57_LbraziliensisMHOMBR75M2904_2019.genes.ext.bed
out_file=${cov_file/.cov/.cnv}

## run depth2cnv.py script
/user/antwerpen/205/vsc20587/scratch/leishmania_susl/bin/lbraz_CNV_depth2cnv.py ${cov_file} ${bed_file} ${out_file}

## run for all cov-files in the results directory
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/lbraz_susl_swga/cnv
# for cov_file in *.cov; do echo $cov_file; gzip $cov_file; done
# for cov_file in *.cov.gz; do sbatch --export=cov_file=${cov_file} /user/antwerpen/205/vsc20587/scratch/leishmania_susl/bin/lbraz_CNV_depth2cnv.slurm; done

## run for all cov files of the WGS data of Senne, but start here
## from the .cov file and do the zipping inside the slurm script to
## speed up the proces
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/lbraz_susl_swga/cnv
# for cov_file in *.cov; do sbatch --export=cov_file=${cov_file} /user/antwerpen/205/vsc20587/scratch/leishmania_susl/bin/lbraz_CNV_depth2cnv.slurm; done
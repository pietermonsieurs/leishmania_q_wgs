# leishmania_q_wgs

## input data
* first run: 36 WGS samples from Ldon
    * project number from GenomeScan: 106214
    * download script: [genomescan_download.slurm](genomescan_download.slurm): file sizes are rather small so can be run on the login node, and no need for sbatch command
* concatenate the different fastq-files: 
    * [inputdata_concatenate_fastq.sh](inputdata_concatenate_fastq.sh): also no sbatch needed, simply run from the command line on the login node

## SNP detection
* BWA: [bwa_parallel.slurm](bwa_parallel.slurm)
* GATK: [GATK_refgenome.sh](GATK_refgenome.sh)
* merge all g.vcf files into one overall g.vcf file and do genotyping?
    * inherit from SomaliLeish samples


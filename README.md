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
* merge all g.vcf files into one overall g.vcf file and do genotyping
    * inherit from SomaliLeish samples
    * [GATK_combine_and_genotype.slurm](GATK_combine_and_genotype.slurm)
* do selection of the relevant genes (known for drug resistance etc.) and do SNP effect evaluation
    * first find the counterparts of the Lbraz genes in Ldon: [drugs_get_orthologs.sh](drugs_get_orthologs.sh)
        * this is based on a prior comparison where reciprokal best blast hits are found between Lbraz and Ldon
    * [GATK_snpeff.sh](GATK_snpeff.sh)
    * visualization of the SNPs: [GATK_visualize_variants.R](GATK_visualize_variants.R)

## somy detection
* run deeptools on all bam files; [somy_deeptools.slurm](somy_deeptools.slurm)
* visualize and check PAT vs non-PAT: [somy_deeptools_plot.R](somy_deeptools_plot.R)

## CNV detection
* CNV detection using the samtools depth approach, and afterwards delineate the genes. Compare the sequencing depth for that gene with the median depth for the corresponding chromosome
* run samtools depth: [cnv_samtools_depth.slurm](cnv_samtools_depth.slurm)
* convert the samtools depth to CNV values:
    * core script: [cnv_depth2cnv.py](cnv_depth2cnv.py)
    * wrapper to run for every bam / depth file: [cnv_depth2cnv_wrapper.slurm](cnv_depth2cnv_wrapper.slurm)
* plot the CNV values and do a statistical test: [cvn_depth2cnv.R](cvn_depth2cnv.R)



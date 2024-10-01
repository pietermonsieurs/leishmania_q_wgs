#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=4
#SBATCH --time=22:30:00
#SBATCH --job-name=bwa_gatk
#SBATCH -A ap_itg_mpu


module load GATK/4.3.0.0-GCCcore-11.3.0-Java-11
module load BCFtools/1.18-GCC-12.3.0
module load Java/17.0.6
# ltopica and laetiopica refs

ref_genome=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/refgenomes/Leishmania_donovani_16Nov2015beta.fa
vcf_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/results/bwa/
gvcf_output="${vcf_dir}/combined.g.vcf.gz"


## combine gvcf
cd $vcf_dir
gatk CombineGVCFs  -R ${ref_genome}  $(for vcf in *.vcf; do  echo "-V ${vcf} "; done) -O ${gvcf_output} 


## create index file and do the genotyping
gatk IndexFeatureFile --input ${vcf_dir}/combined.g.vcf.gz
gatk GenotypeGVCFs -R ${ref_genome} -V ${vcf_dir}/combined.g.vcf.gz -O ${vcf_dir}/combined.vcf.gz


## split up in SNP and indels, and do separate 
## filtering at both levels


cd ${vcf_dir}

##### SNPs #####
## select SNPs and do basic filtering on the SNPs based
## on the GATK recommended settings
gatk SelectVariants -V ${vcf_dir}/combined.vcf.gz  --select-type-to-include SNP  -O ${vcf_dir}/combined_snps.vcf.gz

## add filter information
gatk VariantFiltration -V ${vcf_dir}/combined_snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2"  \
    -filter "QUAL < 200.0" --filter-name "QUAL200" \
    -filter "SOR > 3.0" --filter-name "SOR3"  \
    -filter "FS > 60.0" --filter-name "FS60"  \
    -filter "MQ < 40.0" --filter-name "MQ40"  \
    -O ${vcf_dir}/combined_filter_added_snps.vcf.gz


## select only SNPs that passed the filter 
gatk SelectVariants -R ${ref_genome} -V ${vcf_dir}/combined_filter_added_snps.vcf.gz --exclude-filtered true -O ${vcf_dir}/combined.filtered.snps.vcf.gz
zcat ${vcf_dir}/combined.filtered.snps.vcf.gz | grep -v #" | wc



##### indels #####
gatk SelectVariants -V ${vcf_dir}/combined.vcf.gz -select-type INDEL -O ${vcf_dir}/combined.indels.vcf.gz
zcat ${vcf_dir}/${prefix}.indels.vcf.gz | grep -v #" | wc

## add filter information
gatk VariantFiltration  -V ${vcf_dir}/combined.indels.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 200.0" --filter-name "QUAL200" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -O ${vcf_dir}/combined.filter_added.indels.vcf.gz

## select only indels that passed the filter - 5347 passed filter
gatk SelectVariants -R $ref_genome -V ${vcf_dir}/combined.filter_added.indels.vcf.gz --exclude-filtered true -O ${vcf_dir}/combined.filtered.indels.vcf.gz

zcat ${vcf_dir}/combined.filtered.indels.vcf.gz | grep -v #" | wc


## merge the SNPs and Indels again to produce one output file
java -jar ~/data/software/picard/picard.jar  SortVcf I=${vcf_dir}/combined.filtered.snps.vcf.gz I=${vcf_dir}/combined.filtered.indels.vcf.gz O=${vcf_dir}/combined.filtered.vcf.gz




## do some additional filtering based on the missing genotype
#gt_threshold=.01
#maf_threshold=.05


#bcftools filter -i F_MISSING <= $gt_threshold && MAF >= $maf_threshold' /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_plosntds/genomes_elife.outbreak.BPK72.VL9.SriLanka.one_per_isc.infantum.filtered.vcf.gz > /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_plosntds/genomes_elife.outbreak.BPK72.VL9.SriLanka.one_per_isc.infantum.filtered.missing_gt_${gt_threshold}.maf_${maf_threshold}.vcf

#bcftools filter -i F_MISSING <= $gt_threshold && MAF >= $maf_threshold' /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_plosntds/genomes_elife.outbreak.BPK72.VL9.SriLanka.one_per_isc.infantum.filtered.snps.vcf.gz > /user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_plosntds/genomes_elife.outbreak.BPK72.VL9.SriLanka.one_per_isc.infantum.filtered.snps.missing_gt_${gt_threshold}.maf_${maf_threshold}.vcf

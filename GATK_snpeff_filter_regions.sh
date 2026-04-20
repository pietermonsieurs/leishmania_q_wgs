export PATH=/Users/pmonsieurs/programming/software/bcftools/:/Users/pmonsieurs/programming/software/htslib:$PATH

cd /Users/pmonsieurs/programming/leishmania_q_wgs/results/bwa

## first bgzip the vcf file and run indexing
bgzip combined.filtered.snpeff.vcf
tabix -p vcf combined.filtered.snpeff.vcf.gz

## do renaming and reordering of samples
bcftools reheader -s renaming.tsv combined.filtered.snpeff.vcf.gz -o combined.filtered.snpeff.renamed.vcf.gz
tabix -p vcf combined.filtered.snpeff.renamed.vcf.gz

bcftools view --samples $(bcftools query -l combined.filtered.snpeff.renamed.vcf.gz | sort | paste -sd, -) combined.filtered.snpeff.renamed.vcf.gz -Oz -o combined.filtered.snpeff.renamed.sorted.vcf.gz
tabix -p vcf combined.filtered.snpeff.renamed.sorted.vcf.gz

## create header files
echo -e "CHROM\tPOS\tREF\tALT\tINFO\t$(bcftools query -l combined.filtered.snpeff.renamed.sorted.vcf.gz | tr '\n' '\t')" > combined.filtered.snpeff.Ld31_APQ1_region.500flank.tsv
echo -e "CHROM\tPOS\tREF\tALT\tINFO\t$(bcftools query -l combined.filtered.snpeff.renamed.sorted.vcf.gz | tr '\n' '\t')" > combined.filtered.snpeff.Ld23_Hlocus_region.500flank.tsv


## extract APQ1 region with 500 nt flanking regions
bcftools view -r Ld31:9522-11466 combined.filtered.snpeff.renamed.sorted.vcf.gz \
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO[\t%GT]\n' >> combined.filtered.snpeff.Ld31_APQ1_region.500flank.tsv

## extract MRPA / H-locus region with 500 nt flanking regions
bcftools view -r Ld23:95218-100927 combined.filtered.snpeff.renamed.sorted.vcf.gz \
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO[\t%GT]\n' >> combined.filtered.snpeff.Ld23_Hlocus_region.500flank.tsv


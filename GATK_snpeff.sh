

vcf_dir=/Users/pmonsieurs/programming/leishmania_q_wgs/results/bwa/

cd ${vcf_dir}
java -jar ~/programming/software/snpEff/snpEff.jar Leishmania_donovani_16Nov2015beta ${vcf_dir}/combined.filtered.vcf.gz > ${vcf_dir}/combined.filtered.snpeff.vcf

cd ${vcf_dir}
java -jar ~/programming/software/snpEff/snpEff.jar Leishmania_donovani_16Nov2015beta ${vcf_dir}/combined.filtered.vcf.gz > ${vcf_dir}/combined.filtered.snpeff.vcf
module load BioTools

cd /user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/results/bwa/
vcf_file=combined.filtered.snps.vcf.gz

bcftools query -l ${vcf_file} > samples.txt
bcftools query -f '[%GT\t]\n' ${vcf_file} | awk '
BEGIN {
    while ((getline sample < "samples.txt") > 0) {
        samples[++num_samples] = sample;
    }
    close("samples.txt");
}
{
    for (i=1; i<=NF; i++) {
        if ($i == "0/1" || $i == "1/0") {
            counts[i]++;
        }
    }
    for (i=1; i<=NF; i++) {
        if ($i == "1/1" || $i == "1/1") {
            counts_homo[i]++;
        }
    }
}
END {
    for (i=1; i<=num_samples; i++) {
        print samples[i] "\t" counts[i] "\t" counts_homo[i];
    }
}'
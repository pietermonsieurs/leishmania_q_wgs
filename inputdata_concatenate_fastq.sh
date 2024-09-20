#!/usr/bin/batch

out_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/data/fastq

cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/LeishQ_WGS_allison_20240904

for sample in {001..036}
do
    echo ${sample}
    cat *001-${sample}*R1.fastq.gz > ${out_dir}/${sample}_R1.fastq.gz
    cat *001-${sample}*R2.fastq.gz > ${out_dir}/${sample}_R2.fastq.gz
done



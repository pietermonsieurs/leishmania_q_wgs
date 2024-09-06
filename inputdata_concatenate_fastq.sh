#!/usr/bin/batch

out_dir=/user/antwerpen/205/vsc20587/aitg_data/jcdujardin/DIORAPHTE_genomes_aethiopica_20220523/merged/

cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/DIORAPHTE_genomes_aethiopica_20220523/

for sample in {001..022}
do
    echo ${sample}
    cat *001-${sample}*R1.fastq.gz > ${out_dir}/${sample}_R1.fastq.gz
    cat *001-${sample}*R2.fastq.gz > ${out_dir}/${sample}_R2.fastq.gz
    cat *001-${sample}*R3.fastq.gz > ${out_dir}/${sample}_R3.fastq.gz
done





out_dir=/user/antwerpen/205/vsc20587/aitg_data/jcdujardin/SINGLE_SuSL_20230116/merged
cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/SINGLE_SuSL_20230116/

for sample in {001..030}
do
    echo ${sample}
    cat *001-${sample}*R1.fastq.gz > ${out_dir}/${sample}_R1.fastq.gz
    cat *001-${sample}*R2.fastq.gz > ${out_dir}/${sample}_R2.fastq.gz
done

cd /user/antwerpen/205/vsc20587/scratch/leishmania_susl/data




out_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/swga_2023/
cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/SINGLE_sWGA_20230516/

for sample in {001..040}
for sample in {038..040}
for sample in {001..001}
do
    echo ${sample}
    cat *001-${sample}*R1.fastq.gz > ${out_dir}/${sample}_R1.fastq.gz
    cat *001-${sample}*R2.fastq.gz > ${out_dir}/${sample}_R2.fastq.gz
done



## concatenate (or simply copy) all data from the skinslit / biopsy 
## study

out_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/aethiopica_skinslit/
cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/DIORAPHTE_SuSL_aethiopica_20230814

for sample in {001..008}
do
    echo ${sample}
    cat *105717-001-${sample}*R1.fastq.gz > ${out_dir}/105717-001-${sample}_R1.fastq.gz
    cat *105717-001-${sample}*R2.fastq.gz > ${out_dir}/105717-001-${sample}_R2.fastq.gz
done



## concatenate (or simply copy) all data from the SomaliLeish  
## study

out_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/aethiopica_somalileish/
cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/DIORAPHTE_SuSL_aethiopica_20240416/
project_nr=106070

for sample in {001..008}
do
    echo ${sample}
    cat *${project_nr}-001-${sample}*R1.fastq.gz > ${out_dir}/${project_nr}-001-${sample}_R1.fastq.gz
    cat *${project_nr}-001-${sample}*R2.fastq.gz > ${out_dir}/${project_nr}-001-${sample}_R2.fastq.gz
done



## concatenate (or simply copy) all data from the second batch of  
## the spatialCL study

out_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/aethiopica_skinslit/
cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/DIORAPHTE_SuSL_aethiopica_20240816/
project_nr=106211

for sample in {001..015}
do
    echo ${sample}
    cat *${project_nr}-001-${sample}*R1.fastq.gz > ${out_dir}/${project_nr}-001-${sample}_R1.fastq.gz
    cat *${project_nr}-001-${sample}*R2.fastq.gz > ${out_dir}/${project_nr}-001-${sample}_R2.fastq.gz
done
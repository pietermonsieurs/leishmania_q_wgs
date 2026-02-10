## make symbolic links to the aitg data share

src_dir=/user/antwerpen/205/vsc20587/aitg_data/jcdujardin/LeishQ_RNAseq_20250709
target_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_10X_Q/data/bulk/

for fastq_file in $src_dir/*_R1.fastq.gz; do
    base_name=$(basename "$fastq_file")
    sample_name="${base_name%%_R1.fastq.gz}"
    ln -s "$fastq_file" "${target_dir}/${sample_name}_R1.fastq.gz"
    ln -s "${fastq_file/_R1/_R2}" "${target_dir}/${sample_name}_R2.fastq.gz"
done
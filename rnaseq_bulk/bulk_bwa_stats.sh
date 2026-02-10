cd /user/antwerpen/205/vsc20587/scratch/leishmania_10X_Q/results/bulk/bwa/

output_stats_file=/user/antwerpen/205/vsc20587/scratch/leishmania_10X_Q/results/bulk/bwa/bulk_bwa_htseq_stats_summary.csv
echo -e "Sample,Total_Reads,Mapped_Reads,No_feature_Reads,Ambiguous_Reads,Low_Quality_Reads,Not_Aligned_Reads,Non_Unique_Reads,Feature_reads" > ${output_stats_file}


for flagstat_file in *.flagstat; do
    echo "Processing file: ${flagstat_file}"
    
    ## extract input file
    sample_prefix=${flagstat_file%.flagstat}
    htseq_file=${sample_prefix}.htseq.strand_no.csv
 
    ## extract stats from the flagstat file
    total_reads=$(grep "0 paired in sequencing" ${flagstat_file} | awk '{print $1}')
    mapped_reads=$(grep "mapped (" ${flagstat_file} | awk '{print $1}')
    total_paired_reads

    ## extract from the the htseq ouput file the non-mapping reads (starting
    ## with __, like __no_feature, __ambiguous, etc)
    no_feature_reads=$(grep "^__no_feature" ${htseq_file} | awk '{print $2}')
    ambigous_reads=$(grep "^__ambiguous" ${htseq_file} | awk '{print $2}')
    low_quality_reads=$(grep "^__too_low_aQual" ${htseq_file} | awk '{print $2}')
    not_aligned_reads=$(grep "^__not_aligned" ${htseq_file} | awk '{print $2}')
    non_unique_reads=$(grep "^__alignment_not_unique" ${htseq_file} | awk '{print $2}')

    ## count the sum of the reads mapping to features (i.e. excluding the __no_feature, __ambiguous, etc)
    feature_reads=$(awk '$1 !~ /^__/ { sum += $2 } END { print sum }' ${htseq_file})

    ## write to output file
    echo -e "${sample_prefix},${total_reads},${mapped_reads},${no_feature_reads},${ambigous_reads},${low_quality_reads},${not_aligned_reads},${non_unique_reads},${feature_reads}" >> ${output_stats_file}  

done
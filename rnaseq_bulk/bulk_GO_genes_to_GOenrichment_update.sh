#### bash commands ####

goa_dir=/Users/pmonsieurs/programming/leishmania_10X_Q/results/bulk/bwa_htseq/go/
go_db_dir=/Users/pmonsieurs/programming/data/GO/

cd /Users/pmonsieurs/programming/leishmania_10X_Q/results/bulk/bwa_htseq/go

for genes_file in *.txt; do
# for genes_file in intersect*.txt; do
    ## create output file
    output_file_go=${genes_file/.txt/.GO.csv}
    output_file_kegg=${genes_file/.txt/.KEGG.csv}

    /Users/pmonsieurs/programming/leishmania_10X_Q/bin/bulk_GO_genes_to_GOenrichment.py \
        -g ${goa_dir}/${genes_file} \
        -d ${go_db_dir}/TriTrypDB-48_LdonovaniBPK282A1_GO.gmt \
        -o ${goa_dir}/${output_file_go}

    /Users/pmonsieurs/programming/leishmania_10X_Q/bin/bulk_GO_genes_to_GOenrichment.py \
        -g ${goa_dir}/${genes_file} \
        -d ${go_db_dir}/TriTrypDB-48_LdonovaniBPK282A1.KEGG.gmt \
        -o ${goa_dir}/${output_file_kegg}

done




#### bash commands for the mfuzz clustering output ####

goa_dir=/Users/pmonsieurs/programming/leishmania_10X_Q/results/bulk/bwa_htseq/mfuzz/
go_db_dir=/Users/pmonsieurs/programming/data/GO/


for genes_file in genes*.txt; do
    ## create output file
    output_file_go=${genes_file/.txt/.GO.csv}
    output_file_kegg=${genes_file/.txt/.KEGG.csv}

    /Users/pmonsieurs/programming/leishmania_10X_Q/bin/bulk_GO_genes_to_GOenrichment.py \
        -g ${goa_dir}/${genes_file} \
        -d ${go_db_dir}/TriTrypDB-48_LdonovaniBPK282A1_GO.gmt \
        -o ${goa_dir}/${output_file_go}

    /Users/pmonsieurs/programming/leishmania_10X_Q/bin/bulk_GO_genes_to_GOenrichment.py \
        -g ${goa_dir}/${genes_file} \
        -d ${go_db_dir}/TriTrypDB-48_LdonovaniBPK282A1.KEGG.gmt \
        -o ${goa_dir}/${output_file_kegg}

done

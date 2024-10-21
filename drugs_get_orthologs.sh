
## file names
genes_file=/Users/pmonsieurs/programming/leishmania_susl/data/genes_resistance_Ldon.txt
search_file=/Users/pmonsieurs/programming/leishmania_10X/data/blast/Mapping_TriTrypDB_vs_Sanger.csv


# Loop over the first and second columns of the genes file
while IFS=',' read -r gene_id gene_name; do
    # echo "Searching for gene: $gene_id"
    
    # Use grep to search the second file for each gene_id
    sanger_id=$(grep "${gene_id}.1.1" "$search_file" | cut -f2)
    
    # Print gene_id, sanger_id, and gene_name
    echo "${gene_id},${sanger_id},${gene_name}"
    
done < <(cut -f1,2 -d "," "$genes_file")
#!/usr/bin/env python3

import argparse
import sys

my_debug = 0
sanger_to_tritrypdb_mapping_file = '/Users/pmonsieurs/programming/leishmania_10X/data/blast/Mapping_Sanger_vs_TriTrypDB.csv'

def get_genes_from_file (file):
    print(file)
    gene_fh = open(file, 'r')
    genes = {}
    for line in gene_fh:
        gene = line.split()[0]
        genes[gene] = 1

    return genes.keys()

def get_go_db (file):
    go_fh = open(file, 'r')
    
    go_name = {}
    go_genes = {}

    all_genes_with_a_GO = [""]
    for line in go_fh:
        line = line.rstrip() 
        data = line.split("\t")
        my_debug and print(data[0])
        go_name[data[0]] = data[1]
        go_genes[data[0]] = data[2:len(data)]
        all_genes_with_a_GO = all_genes_with_a_GO + data[2:len(data)]

    # my_debug and print(all_genes_with_a_GO)
    total_number_of_genes_with_a_go = len(set(all_genes_with_a_GO))
    return go_name, go_genes, total_number_of_genes_with_a_go

def gene_convert_to_TriTrypDB (gene_list):
    mapping_fh = open(sanger_to_tritrypdb_mapping_file, 'r')
    mapping_dict = {}
    for line in mapping_fh:
        line = line.rstrip()
        sanger_gene, tritryp_gene = line.split("\t")
        sanger_gene = sanger_gene[:-2] # remove the .1
        tritryp_gene = tritryp_gene[:-4] # remove the .1.1
        mapping_dict[sanger_gene] = tritryp_gene
    
    gene_list_converted = []
    for gene in gene_list:
        if not gene in mapping_dict:
            print(f"gene {gene} not detected in mapping dictionary")
            continue
        gene_converted = mapping_dict[gene]
        gene_list_converted.append(gene_converted)
        my_debug and print(f"gene {gene} -- converted to -- {gene_converted}")

    return gene_list_converted

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description="Start from a list of genes and lookup the GO classes, which can be used as input for GO enrichment stats, to be performed in R")
    parser.add_argument('-g', '--gene_file', required=True, type=str,
                        help="list of Ldon genes for which GO enrichment should be calculated. Gene ID should be in the first column, and file should be tab-separated")
    parser.add_argument('-d', '--db_file', required=True, type=str,
                        help="database containing the GO classes, and a list of genes per GO class. Should be in .gmt format") 
    parser.add_argument('-o', '--out_file', required=True, type=str,
                        help="list of GO classes found back in the gene list, together with all relevant number for calculating a hypergeometric disribution") 

    args = parser.parse_args()
    gene_file = args.gene_file
    db_file = args.db_file
    out_file = args.out_file

    # extract the number of DE genes, and count how many genes are DE. THis
    # value might be needed to calculate the hypergeometric distribution. However,
    # it would be more fair to include only the genes that have a GO assignment. 
    # but in that case, we also have to take as genome size only those genes that 
    # actually have a go assignment. 
    gene_list = get_genes_from_file(gene_file)
    my_debug and print(gene_list)
    my_debug and print(f"number of genes {len(gene_list)}")

    if "TriTryp" in db_file and "Ldonovani" in db_file:
        print("need for conversion of the Sanger GeneIDs to TriTrypDB IDs") 
        gene_list = gene_convert_to_TriTrypDB(gene_list)

    my_debug and print(gene_list)
    number_of_DE_genes = len(gene_list)
    print(f"number of DE genes = {number_of_DE_genes}")

    # find the list of genes per GO, and extract the name of the GO class
    # Addtion
    go_names, go_genes, total_number_of_genes_with_a_go = get_go_db(db_file)
    my_debug and print(go_names)
    my_debug and print(go_genes)

    # genome_size. For now hard coded. could be retrieved based on a gff-file. 
    # should be part of the input parameters?
    genome_size = 8240

    # calculate overlap between both lists, i.e. how many of the DE genes can
    # be found back in a GO class. Write to an ouput file
    out_fh = open(out_file,'w')
    header = f"GO_id;GO_description;nr_found;gene_list;nr_in_GO;nr_sampled;genome_size\n"
    out_fh.write(header)
    for go in go_names:
        my_debug and print(go)
        my_debug and print(set(go_genes[go]))
        my_debug and print(set(gene_list))

        overlap_genes = set(go_genes[go]) & set(gene_list)
        overlap_genes_count = len(overlap_genes)
        my_debug and print(overlap_genes_count)
        if overlap_genes_count > 0:
            my_debug and print(f"GO {go}: go overlap: {overlap_genes}")
            line_out = f"{go};{go_names[go]};{overlap_genes_count};{"|".join(overlap_genes)};{len(go_genes[go])};{number_of_DE_genes};{genome_size}\n"
            print(line_out)
            out_fh.write(line_out)
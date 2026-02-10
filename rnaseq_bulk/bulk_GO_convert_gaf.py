#!/usr/bin/env python3

import sys

input_dir = '/Users/pmonsieurs/programming/leishmania_10X_Q/data/refgenome/'

## two different GAF files available, one raw one, and one curated. Make selection
## in the lines below
gaf_input = f'{input_dir}/TriTrypDB-68_LdonovaniBPK282A1_GO.gaf'
gaf_input = f'{input_dir}/TriTrypDB-68_LdonovaniBPK282A1_Curated_GO.gaf'

gaf_output = gaf_input.replace('BPK282A1', 'BPK282A1v2')
mapping_file = '/Users/pmonsieurs/programming/leishmania_10X_Q/data/refgenome/Mapping_TriTrypDB_vs_Sanger.csv'

open(mapping_file, 'r')
mapping_dict = {}
for line in open(mapping_file, 'r'):
    line = line.rstrip()
    tritryp_gene, sanger_gene = line.split("\t")
    sanger_gene = sanger_gene[:-2] # remove the .1
    tritryp_gene = tritryp_gene[:-4] # remove the .1.1
    print(f"Mapping {sanger_gene} to {tritryp_gene}")
    mapping_dict[tritryp_gene] = sanger_gene
    
# sys.exit()

gaf_fh = open(gaf_input, 'r')
gaf_out_fh = open(gaf_output, 'w')

for line in gaf_fh:
    line = line.rstrip()
    if line.startswith('!'):
        gaf_out_fh.write(line + "\n")
        continue

    data = line.split("\t")
    sanger_gene = data[1]
    sanger_gene = sanger_gene[:-2] # remove the .1

    if not sanger_gene in mapping_dict:
        print(f"gene {sanger_gene} not detected in mapping dictionary")
        continue
    tritryp_gene = mapping_dict[sanger_gene]
    data[1] = tritryp_gene
    new_line = "\t".join(data)
    gaf_out_fh.write(new_line + "\n")

gaf_fh.close()
gaf_out_fh.close()



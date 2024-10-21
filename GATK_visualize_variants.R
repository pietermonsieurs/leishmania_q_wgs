library('vcfR')
library('adegenet')
library('pheatmap')
library('stringr')
library('RColorBrewer')
library('reshape2')
library('ggplot2')
library(openxlsx)
library(gridExtra)
library(viridis)

meta_data_file = '/Users/pmonsieurs/Library/CloudStorage/OneDrive-ITG/leishmania_q_wgs/data/WGS_Samples_DataBase.xlsx'
src_dir = '/Users/pmonsieurs/programming/leishmania_q_wgs/results/bwa/'
setwd(src_dir)


vcf_file = 'combined.filtered.vcf.gz'
vcf_file = 'combined.filtered.snpeff.vcf'





#### 1. overall analysis ####


snp_colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(3)
breaks = c(0,.66,1.33,2.00)

# ---- 1.1 input data ----
## do raw analysis on the genotype as predicted in the .vcf file obtained from
## GATK analsysi
gt <- extract.gt(vcf, element = "GT")
head(gt)


## use vcfR to read data
vcf <- read.vcfR(vcf_file, verbose = TRUE)
queryMETA(vcf)
queryMETA(vcf, element = 'ANN')
head(getFIX(vcf))
vcf@gt[1:5,]

## us genlight object
x <- vcfR2genlight(vcf)
## !! one SNP locus = multiallelic --> only take into account biallelic. So one
## !! locus is omitted in the further analysis for SNP
gt_matrix = t(as.data.frame(x))


## read in meta data. ID_code_short is now a separate column in the meta
## data excel file, and does not need to be derived. Needed as it now contains 
## e.g. 001_Yeti
meta_data = read.xlsx(meta_data_file, startRow = 2, sheet = "Data_base")
# meta_data$ID_code_short = gsub("106214-001-", "", meta_data$ID_code)
meta_data$sample_name = paste0(meta_data$strain, "_", meta_data$PAT, "_", meta_data$Replicate)


## mapping of sample names to biologically interpretable names accoring to meta
## data supplied by Allison
new_col_names = meta_data[match(colnames(gt_matrix), meta_data$ID_code_short),]$sample_name
colnames(gt_matrix) = new_col_names
head(gt_matrix)

## clean up the gt matrix and remove all SNPs where the sum over the different
## samples is 6 (sample 001, 002, 003, 019, 020 and 021). This is BPK026 before
## and after PAT, and contains a very high amount of SNPs compared to the other 
## strains (BPK026 = Yeti, other = ISC core group)
gt_matrix_sub = gt_matrix[- which(rowSums(gt_matrix == 2, na.rm=TRUE) == 6),]

gt_matrix_sub[is.na(gt_matrix_sub)] = 0

## remove the Yeti strains as they contribute with too many SNPs 
## making it impossible to create a heatmap
gt_matrix_sub_select = gt_matrix_sub[, - grep("BPK026", colnames(gt_matrix_sub))]
gt_matrix_sub_select = gt_matrix_sub[, - grep("BPK031", colnames(gt_matrix_sub))]
gt_matrix_sub_select = gt_matrix_sub[, - grep("BPK156", colnames(gt_matrix_sub))]


gt_matrix_sub_select = gt_matrix_sub_select[-which(rowSums(gt_matrix_sub_select == 0) == 42),]
gt_matrix_sub_select = gt_matrix_sub_select[-which(rowSums(gt_matrix_sub_select == 1) == 42),]
gt_matrix_sub_select = gt_matrix_sub_select[-which(rowSums(gt_matrix_sub_select == 2) == 42),]
# gt_matrix_sub_select = gt_matrix_sub_select[-which(rowSums(is.na(gt_matrix_sub)) == 30),]

## add annotation data for heatmap
annotation_data = data.frame(strain = meta_data$strain, PAT = meta_data$PAT)
annotation_data = annotation_data[-grep("BPK026", annotation_data$strain),]
annotation_data = annotation_data[-grep("BPK031", annotation_data$strain),]
annotation_data = annotation_data[-grep("BPK156", annotation_data$strain),]

rownames(annotation_data) = colnames(gt_matrix_sub_select)



pheatmap(gt_matrix_sub_select,
         fontsize_row = 2,
         show_rownames = FALSE,
         cluster_cols = TRUE,
         color = snp_colors,
         breaks = breaks, 
         treeheight_row = FALSE,
         annotation_col = annotation_data)


# ---- 1.2 focus on SNPs in coding regions ----

## only select those position where the SNP has in impact, which means that it 
## should occur in the coding region, and the SNP should have an impact (i.e. 
## a nonsense or missense mutation). The impact of those SNPs and indels can be
## assessed using the SnpEff approach, 

## read in the info field from the vcf file, which contains the information of 
## SNPeff
info_data = INFO2df(vcf)
ann_data = str_split_fixed(info_data$ANN, "\\|", n=10)
head(ann_data)
table(ann_data[,2])
table(ann_data[,5])


## select only those SNP with a substantial impact (e.g. ("missense_variant", 
## "stop_gained") or remove the one with a non-harmful deletion (e.g. 
## "downstream_gene_variant", "synonymous_variant", "upstream_gene_variant")
variants_omit = c("downstream_gene_variant", "synonymous_variant", "upstream_gene_variant", "intergenic_region")
# mutation_impact = ann_data[,2] %in% variants_omit
# vcf_impact = vcf[-mutation_impact,]
vcf_impact = vcf[! ann_data[,2] %in% variants_omit, ]


# 
# if (type == 'snp') {
#   variants_select = c("missense_variant", "stop_gained")
#   mutation_impact = ann_data[,2] %in% variants_select
# }else if (type == 'indel') {
#   variants_select = c("conservative_inframe_deletion", "frameshift_variant")
#   mutation_impact = ann_data[,2] %in% variants_select
# }
#   
# vcf_impact = vcf[mutation_impact,]


info_data_genename = INFO2df(vcf_impact)
ann_data_genename = str_split_fixed(info_data_genename$ANN, "\\|", n=10)
head(ann_data_genename)
table(ann_data_genename[,2])
table(ann_data_genename[,5])

dim(ann_data_genename)
dim(vcf_impact@fix)
head(vcf_impact@fix)

## recreate the genename in the annotation data frame by concatenating the
## chromosome with the position
rownames(ann_data_genename) = paste0(vcf_impact@fix[,1], "_", vcf_impact@fix[,2])
head(ann_data_genename[,1:5])




## create heatmap with impact SNPs
x_impact <- vcfR2genlight(vcf_impact)
gt_impact_matrix = t(as.data.frame(x_impact))

## mapping of sample names to biologically interpretable names accoring to meta
## data supplied by Allison

colnames(gt_impact_matrix) = new_col_names
head(gt_matrix)

## clean up the gt matrix and remove all SNPs where the sum over the different
## samples is 6 (sample 001, 002, 003, 019, 020 and 021). This is BPK026 before
## and after PAT, and contains a very high amount of SNPs compared to the other 
## strains (BPK026 = Yeti, other = ISC core group)
# gt_matrix_sub = gt_impact_matrix[- which(rowSums(gt_impact_matrix == 2, na.rm=TRUE) == 6),]
gt_matrix_sub = gt_impact_matrix

gt_matrix_sub[is.na(gt_matrix_sub)] = 0

## remove samples 001,002 etc.
# gt_matrix_sub_select = gt_matrix_sub[, - which(colnames(gt_matrix_sub) %in% c('001', '002', '003', '019', '020', '021'))]
gt_matrix_sub_select = gt_matrix_sub[, - grep("BPK026", colnames(gt_matrix_sub))]
gt_matrix_sub_select = gt_matrix_sub[, - grep("BPK031", colnames(gt_matrix_sub))]
gt_matrix_sub_select = gt_matrix_sub[, - grep("BPK156", colnames(gt_matrix_sub))]

gt_matrix_sub_select = gt_matrix_sub_select[-which(rowSums(gt_matrix_sub_select == 0) == 42),]
gt_matrix_sub_select = gt_matrix_sub_select[-which(rowSums(gt_matrix_sub_select == 1) == 42),]
# gt_matrix_sub_select = gt_matrix_sub_select[-which(rowSums(gt_matrix_sub_select == 2) == 30),]
# gt_matrix_sub_select = gt_matrix_sub_select[-which(rowSums(is.na(gt_matrix_sub)) == 30),]

## add annotation data for heatmap
annotation_data = data.frame(strain = meta_data$strain, PAT = meta_data$PAT)
annotation_data = annotation_data[-grep("BPK026", annotation_data$strain),]
annotation_data = annotation_data[-grep("BPK031", annotation_data$strain),]
annotation_data = annotation_data[-grep("BPK156", annotation_data$strain),]

rownames(annotation_data) = colnames(gt_matrix_sub_select)



pheatmap(gt_matrix_sub_select,
         fontsize_row = 2,
         show_rownames = FALSE,
         cluster_cols = TRUE,
         color = snp_colors,
         breaks = breaks, 
         treeheight_row = FALSE,
         annotation_col = annotation_data)


# ---- 1.3 remove overall heterozyogtes ----
## remove all the positions that are heterozygous for all positions
# gt_matrix_filt = gt_impact_matrix[rowSums(gt_impact_matrix == 1) < 6,]
# pheatmap(gt_matrix_filt,
#          fontsize_row = 10,
#          cluster_cols = FALSE,
#          treeheight_row = FALSE)
# 
# ## remove all the positions that contain at least one position with a 
# ## homozygote alternative allele
# gt_matrix_filt2 = gt_matrix[rowSums(gt_matrix == 2) > 0,]
# pheatmap(gt_matrix_filt2,
#          fontsize_row = 16,
#          fontsize_col = 16,
#          cluster_cols = FALSE,
#          treeheight_row = FALSE)



# ---- 1.4 analysis per strain ----

strains = unique(meta_data$strain)
p_heatmaps = list()
p_heatmaps_impact = list()
p_heatmaps_impact_homo = list()

meta_data_sorted = meta_data[match(colnames(gt_matrix), meta_data$sample_name),]
annotation_data = data.frame(strain = meta_data_sorted$strain, PAT = meta_data_sorted$PAT)
# annotation_data = annotation_data[-grep("BPK026", annotation_data$strain),]
rownames(annotation_data) = colnames(gt_matrix)


for (strain in strains) {
  head(gt_matrix)
  
  ## analysis for all SNPs combined
  gt_matrix_strain = gt_matrix[,grep(strain, colnames(gt_matrix))]
  # gt_matrix_strain = gt_matrix_strain[- which(rowSums(gt_matrix_strain == 2, na.rm=TRUE) == 6),]
  gt_matrix_strain[is.na(gt_matrix_strain)] = 0
  gt_matrix_strain_cleaned <- gt_matrix_strain[!apply(gt_matrix_strain, 1, function(row) length(unique(row)) == 1), ]
  
    
  p = pheatmap(gt_matrix_strain_cleaned,
           fontsize_row = 2,
           show_rownames = FALSE,
           cluster_cols = FALSE,
           color = snp_colors,
           breaks = breaks, 
           treeheight_row = FALSE,
           treeheight_col = FALSE,
           annotation_col = annotation_data)
  p_heatmaps[[strain]] = p
  
  ## for SNPs with an impact on the CDS
  gt_impact_matrix_strain = gt_impact_matrix[,grep(strain, colnames(gt_impact_matrix))]
  gt_impact_matrix_strain[is.na(gt_impact_matrix_strain)] = 0
  gt_impact_matrix_strain_cleaned <- gt_impact_matrix_strain[!apply(gt_impact_matrix_strain, 1, function(row) length(unique(row)) == 1), ]
  
  
  p_impact = pheatmap(gt_impact_matrix_strain_cleaned,
               fontsize_row = 2,
               show_rownames = FALSE,
               cluster_cols = FALSE,
               color = snp_colors,
               breaks = breaks, 
               treeheight_row = FALSE,
               treeheight_col = FALSE,
               annotation_col = annotation_data)
  p_heatmaps_impact[[strain]] = p_impact
  
  
  
  gt_impact_matrix_strain_cleaned_homo = gt_impact_matrix_strain_cleaned[rowSums(gt_impact_matrix_strain_cleaned == 2) > 0,]
  
  p_impact_homo = pheatmap(gt_impact_matrix_strain_cleaned_homo,
                      fontsize_row = 2,
                      show_rownames = FALSE,
                      cluster_cols = FALSE,
                      color = snp_colors,
                      breaks = breaks, 
                      treeheight_row = FALSE,
                      treeheight_col = FALSE,
                      annotation_col = annotation_data)
  p_heatmaps_impact_homo[[strain]] = p_impact_homo
  
  
}


grid.arrange(p_heatmaps[["BPK026"]]$gtable, 
             p_heatmaps[["BPK031"]]$gtable, 
             p_heatmaps[["BPK156"]]$gtable, 
             p_heatmaps[["BPK275"]]$gtable, 
             p_heatmaps[["BPK080"]]$gtable, 
             p_heatmaps[["BPK085"]]$gtable, 
             p_heatmaps[["BPK282"]]$gtable, 
             p_heatmaps[["BPK294"]]$gtable, 
             ncol = 4)

grid.arrange(p_heatmaps_impact[["BPK026"]]$gtable, 
             p_heatmaps_impact[["BPK031"]]$gtable, 
             p_heatmaps_impact[["BPK156"]]$gtable, 
             p_heatmaps_impact[["BPK275"]]$gtable, 
             p_heatmaps_impact[["BPK080"]]$gtable, 
             p_heatmaps_impact[["BPK085"]]$gtable, 
             p_heatmaps_impact[["BPK282"]]$gtable, 
             p_heatmaps_impact[["BPK294"]]$gtable, 
             ncol = 4)

grid.arrange(p_heatmaps_impact_homo[["BPK026"]]$gtable,
             p_heatmaps_impact_homo[["BPK031"]]$gtable,
             p_heatmaps_impact_homo[["BPK156"]]$gtable,
             p_heatmaps_impact_homo[["BPK275"]]$gtable, 
             p_heatmaps_impact_homo[["BPK080"]]$gtable, 
             p_heatmaps_impact_homo[["BPK085"]]$gtable, 
             p_heatmaps_impact_homo[["BPK282"]]$gtable, 
             p_heatmaps_impact_homo[["BPK294"]]$gtable, 
             ncol = 4)

p_heatmaps[['BPK026']]


## overlap of SNPs between different strains

## intialize empty dataframe
p_overlap = list()

for (strain in strains) {

  ## for SNPs with an impact on the CDS
  gt_impact_matrix_strain = gt_impact_matrix[,grep(strain, colnames(gt_impact_matrix))]
  gt_impact_matrix_strain[is.na(gt_impact_matrix_strain)] = 0
  gt_impact_matrix_strain_cleaned <- gt_impact_matrix_strain[!apply(gt_impact_matrix_strain, 1, function(row) length(unique(row)) == 1), ]
  
  sample_names = colnames(gt_impact_matrix_strain_cleaned)
  snp_overlap = data.frame(matrix(0, nrow = length(sample_names), ncol = length(sample_names)))
  colnames(snp_overlap) = sample_names
  rownames(snp_overlap) = sample_names
  
  for (sample_nr in 1:length(sample_names)) {
    for (sample_nr2 in 1:length(sample_names)) {
      sample1 = sample_names[sample_nr]
      sample2 = sample_names[sample_nr2]
      overlap = sum(gt_impact_matrix_strain_cleaned[,sample1] > 0 & gt_impact_matrix_strain_cleaned[,sample2] > 0)
      snp_overlap[sample1, sample2] = overlap
    }
  }

  
  p = pheatmap(snp_overlap,
                      color = rev(inferno(20)),
                      fontsize_row = 2,
                      show_rownames = FALSE,
                      cluster_cols = FALSE,
                      cluster_rows = FALSE,
                      treeheight_row = FALSE,
                      treeheight_col = FALSE,
                      display_numbers = TRUE,
                      number_format = "%.0f",
                      fontsize_number = 18,
                      number_color = "white",
                      main = strain)
             
                    
                      
  p_overlap[[strain]] = p
  
  
}

grid.arrange(p_overlap[["BPK026"]]$gtable,
             p_overlap[["BPK031"]]$gtable,
             p_overlap[["BPK156"]]$gtable,
             p_overlap[["BPK275"]]$gtable, 
             p_overlap[["BPK080"]]$gtable, 
             p_overlap[["BPK085"]]$gtable, 
             p_overlap[["BPK282"]]$gtable, 
             p_overlap[["BPK294"]]$gtable, 
             ncol = 4)


## combine all Yeti strains



#### 2. summarise SNPs per gene ####

# ---- 2.1 check gene information ---- #
head(gt_impact_matrix)
head(ann_data_genename)
dim(gt_impact_matrix)
dim(ann_data_genename)
rownames(ann_data_genename)

## make a copy of gt_impact_matrix
gt_impact_matrix_copy = gt_impact_matrix
gt_impact_matrix_copy <- as.data.frame(gt_impact_matrix_copy)
gene_names = ann_data_genename[match(rownames(gt_impact_matrix_copy), rownames(ann_data_genename)),5]
gt_impact_matrix_copy$gene = unlist(gene_names, use.names = FALSE)


## summarize all SNPs per gene and select the maximum value (which can 
## be 0, 1 or 2)
library(dplyr)
gt_impact_matrix_gene <- gt_impact_matrix_copy %>%
  group_by(gene) %>%
  summarise(across(everything(), max, na.rm = TRUE))


## replace all values with a value of 2 by 1, as the first plan is to 
## not make a distinction between homozygote and heterozygote SNPs
gt_impact_matrix_gene[gt_impact_matrix_gene == 2] <- 1

## add the gene as rowname and remove it from the data frame
gt_impact_matrix_gene = as.data.frame(gt_impact_matrix_gene)
rownames(gt_impact_matrix_gene) = gt_impact_matrix_gene$gene
gt_impact_matrix_gene = gt_impact_matrix_gene[,-1]
head(gt_impact_matrix_gene)
dim(gt_impact_matrix_gene)

## do the distinction between Yeti strains and core strains.  
# rowSums(gt_impact_matrix_gene)
colnames(gt_impact_matrix_gene)
yeti_columns = grep("BPK026|BPK031|BPK156", colnames(gt_impact_matrix_gene))

## yeti-analysis: select only the Yeti strains or *not* on Yeti strains
gt_impact_matrix_gene_yeti = gt_impact_matrix_gene[, yeti_columns]
gt_impact_matrix_gene_yeti = gt_impact_matrix_gene[, - yeti_columns]


gt_impact_matrix_gene_yeti[gt_impact_matrix_gene_yeti == -Inf] <- NA
dim(gt_impact_matrix_gene_yeti)
gt_impact_matrix_gene_yeti <- gt_impact_matrix_gene_yeti[!apply(gt_impact_matrix_gene_yeti, 1, function(row) length(unique(row)) == 1), ]
dim(gt_impact_matrix_gene_yeti)

# Replace NA or Inf with 0
gt_impact_matrix_gene_yeti[is.na(gt_impact_matrix_gene_yeti)] <- 0
#gt_impact_matrix_gene_yeti[!is.finite(gt_impact_matrix_gene_yeti)] <- 0  # Handle Inf/-Inf
gt_impact_matrix_gene_yeti[!apply(gt_impact_matrix_gene_yeti, 2, is.finite)] <- 0  # Handle Inf/-Inf


annotation_data = as.data.frame(meta_data[match(colnames(gt_impact_matrix_gene_yeti), meta_data$sample_name),]$PAT)
rownames(annotation_data) = colnames(gt_impact_matrix_gene_yeti)
colnames(annotation_data) = c('PAT')

library(pheatmap)
pheatmap(gt_impact_matrix_gene_yeti,
         cluster_cols = TRUE,
         cluster_rows = TRUE, 
         fontsize_row = 5, 
         annotation_col = annotation_data)


## remove the strain specific SNPs as they will not be contributing to 
## the drug resistance
strains = unique(meta_data$strain)
rows_all_remove = c()
for (strain in strains) {
  gt_impact_matrix_gene_yeti_strain = gt_impact_matrix_gene_yeti[, grep(strain, colnames(gt_impact_matrix_gene_yeti))]
  rows_to_remove = which(apply(gt_impact_matrix_gene_yeti_strain, 1, function(row) length(unique(row)) == 1))
  rows_all_remove = c(rows_all_remove, rows_to_remove)
}

gt_impact_matrix_gene_yeti_clean = gt_impact_matrix_gene_yeti[-rows_all_remove,]
pheatmap(gt_impact_matrix_gene_yeti_clean,
         cluster_cols = TRUE,
         cluster_rows = TRUE, 
         fontsize_row = 5, 
         annotation_col = annotation_data)





#### 3. drug resistance regions ####

## make a copy of gt_impact_matrix
gt_impact_matrix_copy = gt_impact_matrix
gt_impact_matrix_copy <- as.data.frame(gt_impact_matrix_copy)
gene_names = ann_data_genename[match(rownames(gt_impact_matrix_copy), rownames(ann_data_genename)),5]
gt_impact_matrix_copy$gene = unlist(gene_names, use.names = FALSE)

## read in the excel file containing the drug resistance regions and extract
## them from the bigger file 
dr_genes_file = '/Users/pmonsieurs/programming/leishmania_q_wgs/data/drug_genes_Ldon.xlsx'
drug_genes = read.xlsx(dr_genes_file)
drug_genes$gene_id_v2 = gsub(".1$", "", drug_genes$gene_id_v2)
head(drug_genes)

## select only the drug resistance genes from the gt_matrix
gt_impact_matrix_dr = gt_impact_matrix_copy
gt_impact_matrix_dr = gt_impact_matrix_dr[gt_impact_matrix_dr$gene %in% drug_genes$gene_id_v2,]
gt_impact_matrix_dr$drug_gene = drug_genes[match(gt_impact_matrix_dr$gene, drug_genes$gene_id_v2),]$resistance
                                      
## get the snpEff information 
snp_eff_data = ann_data_genename[match(rownames(gt_impact_matrix_dr), rownames(ann_data_genename)),c(2,3)]
colnames(snp_eff_data) = c('snp_variant', 'snp_impact')
if (dim(snp_eff_data)[1] == dim(gt_impact_matrix_dr)[1]) {
  gt_impact_matrix_dr = cbind.data.frame(gt_impact_matrix_dr, snp_eff_data)
}

gt_impact_matrix_dr$snp_variant = gsub("_variant", "", gt_impact_matrix_dr$snp_variant)
gt_impact_matrix_dr$rowname = paste0(gt_impact_matrix_dr$drug_gene, '|', 
                                     rownames(gt_impact_matrix_dr), '|',
                                     gt_impact_matrix_dr$snp_variant)
rownames(gt_impact_matrix_dr) = gt_impact_matrix_dr$rowname


pheatmap(gt_impact_matrix_dr[,1:48],
         cluster_rows=FALSE,
         cluster_cols=TRUE,
         treeheight_col = 0)





####### BELOW NOT (YET?) USED FOR Q WGS ######




#### 2. Allele frequencies ####

# ---- 2.1 extract allele depths ----

## use the extract.gt to extract allele depths
ad <- extract.gt(vcf, element = "AD")
ad1 <- masplit(ad, record = 1)
ad2 <- masplit(ad, record = 2)
overall_depth = ad1 + ad2

overall_depth_plotdata = melt(overall_depth)
overall_depth_plotdata = melt(ad1)
overall_depth_plotdata = melt(ad2)

colnames(overall_depth_plotdata) = c('snp', 'sample', 'depth')
head(overall_depth_plotdata)

p = ggplot(overall_depth_plotdata, aes(depth, colour=sample))
p = p + geom_density()
p = p + theme_bw()
p


# ---- 2.2 extract allele depths for impact SNPs ----
ad <- extract.gt(vcf_impact, element = "AD")
ad1 <- masplit(ad, record = 1, sort=0) ## sort is needed, otherwise the highest depth is put first
ad2 <- masplit(ad, record = 2, sort=0)
overall_depth = ad1 + ad2
af = ad2/(ad1+ad2)
colnames(af) = c('DMSO', 'C', 'D1', 'D2', 'D3', 'D4')



# ---- 2.3 create heatmap based on allele frequency ----
pheatmap(af,
         fontsize_row = 10,
         cluster_cols = FALSE,
         treeheight_row = FALSE)

af_filtered = af[rowSums(gt_impact_matrix == 1) < 6,]
pheatmap(af_filtered,
         fontsize_row = 10,
         cluster_cols = FALSE,
         treeheight_row = FALSE)

# ---- 2.4 get overall depth matrix ---- 
depth_colors = colorRampPalette(brewer.pal(n = 7, name ="Reds"))(100)
colnames(overall_depth) = c('DMSO', 'C', 'D1', 'D2', 'D3', 'D4')

pheatmap(overall_depth,
         color = depth_colors,
         fontsize_row = 10,
         cluster_cols = FALSE,
         treeheight_row = FALSE)

breaks = seq(0,100,1)
depth_colors = colorRampPalette(brewer.pal(n = 7, name ="Reds"))(length(breaks) -1)


overall_depth_filt = overall_depth[rowSums(gt_impact_matrix == 1) < 6,]
pheatmap(overall_depth_filt,
         color = depth_colors,
         breaks = breaks,
         fontsize_row = 10,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         treeheight_row = FALSE)

## make same heatmaps as above for genotypes and allele frequencies but 
## do not do clustering based on rows, but keep ordering based on chromosomes
gt_matrix_filt = gt_impact_matrix[rowSums(gt_impact_matrix == 1) < 6,]
pheatmap(gt_matrix_filt,
         fontsize_row = 10,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         treeheight_row = FALSE)

af_filtered = af[rowSums(gt_impact_matrix == 1) < 6,]
pheatmap(af_filtered,
         fontsize_row = 10,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         treeheight_row = FALSE)

## map informaition from word file JC
info_data2 = INFO2df(vcf_impact)
ann_data2 = str_split_fixed(info_data2$ANN, "\\|", n=10)
ann_data2 = ann_data2[rowSums(gt_impact_matrix == 1) < 6,]
head(ann_data2)
ann_data2[,4]

## create pheatmap with product name added to the SNP position
if (type == 'indel') {af_filtered = af}
row_labels = paste0(rownames(af_filtered), " || ", ann_data2[,4])
pheatmap(af_filtered,
         fontsize_row = 10,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         treeheight_row = FALSE,
         labels_row = row_labels)

pheatmap(af_filtered,
         fontsize_row = 10,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         treeheight_row = FALSE,
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 10,
         number_color = "black",
         labels_row = row_labels)

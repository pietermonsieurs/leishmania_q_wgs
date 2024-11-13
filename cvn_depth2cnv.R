library(ggplot2)
library(openxlsx)
library(pheatmap)
library(gridExtra)


src_dir = '/Users/pmonsieurs/programming/leishmania_q_wgs/results/cnv/'
cnv_files = list.files(src_dir, pattern = ".cnv.csv")
meta_data_file = '/Users/pmonsieurs/Library/CloudStorage/OneDrive-ITG/leishmania_q_wgs/data/WGS_Samples_DataBase.xlsx'


## read in the meta data and create an additional column name
meta_data = read.xlsx(meta_data_file, startRow = 2, sheet = "Data_base")
meta_data$sample_name = paste0(meta_data$strain, "_", meta_data$PAT, "_", meta_data$Replicate)
head(meta_data)
meta_data$ID_code_short


first = 1
for (cnv_file in cnv_files) {
  cnv_new_data = read.csv(paste0(src_dir, cnv_file))
  
  ## concatenate the new data to the existing data frame
  if (first == 1) {
    cnv_data = cnv_new_data
    first = 0
  }else{
    cnv_data = cbind.data.frame(cnv_data, cnv_new_data$ratio)
  }
  
  ## create the column names
  sample = gsub(".cnv.csv", "", cnv_file)
  sample_name = meta_data[meta_data$ID_code_short == sample,]$sample_name
  colnames(cnv_data)[length(colnames(cnv_data))] = sample_name
  
}


dim(cnv_data)
head(cnv_data)


## filter out those CNV values that are higher than 1.5
cutoff_cnv = 1.5
cnv_data_up = cnv_data[rowSums(cnv_data[,5:ncol(cnv_data)] > cutoff_cnv) > 0,]
rownames(cnv_data_up) = cnv_data_up$gene_id

## create heatmap with all strains combined
breaks = seq(0,5,0.05)
pheatmap(t(cnv_data_up[,5:ncol(cnv_data_up)]),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         treeheight_col = 0,
         breaks = breaks)

## create heatmap per strain
strains = unique(meta_data$strain)
p_heatmaps = list()
p_heatmaps_significant = list()
relevant_genes = c()

for (strain in strains) {
  col_indices =  grep(strain, colnames(cnv_data))
  col_indices = c(1:4, col_indices)
  cnv_data_strain = cnv_data[,col_indices]
  
  cnv_data_strain = cnv_data_strain[rowSums(cnv_data_strain[,5:ncol(cnv_data_strain)] > cutoff_cnv) > 0,]
  
  breaks = seq(0,5,0.05)
  p = pheatmap(t(cnv_data_strain[,5:ncol(cnv_data_strain)]),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_colnames = FALSE,
           treeheight_col = 0,
           treeheight_row = 0,
           breaks = breaks,
           fontsize_row = 8)
  
  p_heatmaps[[strain]] = p
  
  
  ## do t-test per gene
  p_values <- numeric(nrow(cnv_data_strain))
  for (i in 1:nrow(cnv_data_strain)) {
    # Extract values for columns 5 to 8 and 9 to 11 for the current row
    group_noPAT <- as.numeric(cnv_data_strain[i, 5:7])
    group_PAT <- as.numeric(cnv_data_strain[i, 8:10])
    
    # Perform the t-test and store the p-value
    test_result <- t.test(group_noPAT, group_PAT)
    p_values[i] <- test_result$p.value
  }
  
  ## add the p-values as a new column to the cnv_data dataframe, and filter
  ## based on this p-value
  cnv_data_strain$p_value <- p_values
  cnv_data_strain_sign <- cnv_data_strain[cnv_data_strain$p_value < 0.05, ]
  rownames(cnv_data_strain_sign) = cnv_data_strain_sign$gene_id
  
  ## make heatmap of only those CVN's that are significant
  breaks = seq(0,5,0.05)
  p = pheatmap(t(cnv_data_strain_sign[,5:10]),
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_colnames = TRUE,
               treeheight_col = 0,
               treeheight_row = 0,
               breaks = breaks,
               fontsize_row = 8)
  
  p_heatmaps_significant[[strain]] = p
  
  ## add to the relevant genes, that might be interesting to make a separate
  ## heatmap 
  relevant_genes = c(relevant_genes, cnv_data_strain_sign$gene_id)
  
}

grid.arrange(p_heatmaps[["BPK026"]]$gtable, 
             p_heatmaps[["BPK031"]]$gtable, 
             p_heatmaps[["BPK156"]]$gtable, 
             p_heatmaps[["BPK275"]]$gtable, 
             p_heatmaps[["BPK080"]]$gtable, 
             p_heatmaps[["BPK085"]]$gtable, 
             p_heatmaps[["BPK282"]]$gtable, 
             p_heatmaps[["BPK294"]]$gtable, 
             ncol = 3)


grid.arrange(p_heatmaps_significant[["BPK026"]]$gtable, 
             p_heatmaps_significant[["BPK031"]]$gtable, 
             p_heatmaps_significant[["BPK156"]]$gtable, 
             p_heatmaps_significant[["BPK275"]]$gtable, 
             p_heatmaps_significant[["BPK080"]]$gtable, 
             p_heatmaps_significant[["BPK085"]]$gtable, 
             p_heatmaps_significant[["BPK282"]]$gtable, 
             p_heatmaps_significant[["BPK294"]]$gtable, 
             ncol = 3)


## make heatmap with only those genes that were significant somewhere in the 
## statistical test of the different strains
length(relevant_genes)
relevant_genes_unique = unique(relevant_genes)
cnv_data_relevant = cnv_data[cnv_data$gene_id %in% relevant_genes_unique, ]
rownames(cnv_data_relevant) = cnv_data_relevant$gene_id

pheatmap(t(cnv_data_relevant[,5:ncol(cnv_data_relevant)]),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = TRUE,
         treeheight_col = 0,
         treeheight_row = 0,
         breaks = breaks)

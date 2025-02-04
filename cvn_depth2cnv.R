library(ggplot2)
library(openxlsx)
library(pheatmap)
library(gridExtra)
library(reshape2)
library(stringr)
library(patchwork)

src_dir = '/Users/pmonsieurs/programming/leishmania_q_wgs/results/cnv/'
cnv_files = list.files(src_dir, pattern = ".cnv.csv")
meta_data_file = '/Users/pmonsieurs/Library/CloudStorage/OneDrive-ITG/leishmania_q_wgs/data/WGS_Samples_DataBase.xlsx'
out_dir = '/Users/pmonsieurs/programming/leishmania_q_wgs/docs/figures/draft_nov2024/'


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

## same but only for Yeti
col_index = c(grep("BPK026", colnames(cnv_data_relevant)),
              grep("BPK031", colnames(cnv_data_relevant)),
              grep("BPK156", colnames(cnv_data_relevant)))
cnv_data_yeti = cnv_data[, c(1:4,col_index)]
cnv_data_yeti_up = cnv_data_yeti[rowSums(cnv_data_yeti[,5:ncol(cnv_data_yeti)] > cutoff_cnv) > 0,]
cnv_data_yeti_down = cnv_data_yeti[rowSums(cnv_data_yeti[,5:ncol(cnv_data_yeti)] < 1/cutoff_cnv) > 0,]
cnv_data_yeti_diff = rbind.data.frame(cnv_data_yeti_up, cnv_data_yeti_down)
dim(cnv_data_yeti_diff)

breaks = seq(0,5,0.05)
pheatmap(t(cnv_data_yeti_diff[,5:ncol(cnv_data_yeti_diff)]),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         treeheight_col = 0,
         breaks = breaks)


## create heatmap per strain
strains = unique(meta_data$strain)
p_heatmaps = list()
p_heatmaps_significant = list()
p_lineplots = list()
relevant_genes = c()
relevant_genes_core = c()
relevant_genes_yeti = c()
yeti_strains = c('BPK026', 'BPK031', 'BPK156')

for (strain in strains) {
  col_indices =  grep(strain, colnames(cnv_data))
  col_indices = c(1:4, col_indices)
  cnv_data_strain = cnv_data[,col_indices]
  
  # cnv_data_strain = cnv_data_strain[rowSums(cnv_data_strain[,5:ncol(cnv_data_strain)] > cutoff_cnv) > 0,]
  cnv_data_strain_up = cnv_data_strain[rowSums(cnv_data_strain[,5:ncol(cnv_data_strain)] > cutoff_cnv) > 0,]
  cnv_data_strain_down = cnv_data_strain[rowSums(cnv_data_strain[,5:ncol(cnv_data_strain)] < 1/cutoff_cnv) > 0,]
  cnv_data_strain = rbind.data.frame(cnv_data_strain_up, cnv_data_strain_down)
  # cnv_data_strain <- cnv_data_strain %>% distinct()
  
    
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
  p_values[is.na(p_values)] = 1
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
  
  
  ## make line plot 
  plot_data_sign = melt(cnv_data_strain_sign[,c(1,3,5:(ncol(cnv_data_strain_sign)-1))])
  colnames(plot_data_sign) = c('gene_id', 'chrom', 'sample', 'cnv')
  plot_data_sign$sample = as.character(plot_data_sign$sample)
  
  split_col <- sapply(strsplit(plot_data_sign$sample, "_"), function(x) {
    first <- x[1]                     # First element
    middle <- paste(x[2:(length(x) - 1)], collapse = "_")  # Middle elements combined
    last <- x[length(x)]              # Last element
    c(first, middle, last)
  })
  
  split_df = as.data.frame(t(split_col))
  colnames(split_df) = c('strain', 'condition', 'replicate')
  plot_data_sign = cbind(plot_data_sign, split_df)
  
  
  p = ggplot(plot_data_sign, aes(x = gene_id, y = cnv, group = strain)) +
    geom_point(aes(color = condition), size = 3, alpha=0.50) +  # Replicate points colored by Condition
    labs(x = "gene",
         y = "CNV value",
         color = "Condition",
         shape = "Replicate") +
    theme_minimal() + 
    theme(
      axis.text.y = element_text(size = 10)  # Adjust font size of y-axis (row) text
    ) + 
    coord_flip() + 
    facet_wrap(~ strain, scales = c("free_x"), ncol=4)
  
  p_lineplots[[strain]] = p
  
  ## add to the relevant genes, that might be interesting to make a separate
  ## heatmap 
  relevant_genes = c(relevant_genes, cnv_data_strain_sign$gene_id)
  
  if (strain %in% yeti_strains) {
    relevant_genes_yeti = c(relevant_genes_yeti, cnv_data_strain_sign$gene_id)    
  }else{
    relevant_genes_core = c(relevant_genes_core, cnv_data_strain_sign$gene_id)
  }
  
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
grid.arrange(p_lineplots[["BPK026"]], 
             p_lineplots[["BPK031"]], 
             p_lineplots[["BPK156"]], 
             p_lineplots[["BPK275"]], 
             p_lineplots[["BPK080"]], 
             p_lineplots[["BPK085"]], 
             p_lineplots[["BPK282"]], 
             p_lineplots[["BPK294"]], 
             ncol = 3)

## make heatmap with only those genes that were significant somewhere in the 
## statistical test of the different strains

## for ALL genes
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



## for CORE genes

## grep the significant genes
length(relevant_genes_core)
relevant_genes_unique = unique(relevant_genes_core)
cnv_data_relevant = cnv_data[cnv_data$gene_id %in% relevant_genes_unique, ]
rownames(cnv_data_relevant) = cnv_data_relevant$gene_id

## select core strains
col_index = c(grep("BPK026", colnames(cnv_data_relevant)),
              grep("BPK031", colnames(cnv_data_relevant)),
              grep("BPK156", colnames(cnv_data_relevant)))
cnv_data_relevant = cnv_data_relevant[, -col_index]

png_cnv_core = pheatmap(t(cnv_data_relevant[,5:ncol(cnv_data_relevant)]),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = TRUE,
         treeheight_col = 0,
         treeheight_row = 0,
         breaks = breaks)
png_cnv_core

png_cnv_core_file = paste(out_dir, 'fig_04B_CNV_corestrains.png')
ggsave(filename = png_cnv_core_file,
       plot = png_cnv_core,
       width = 16,
       height = 9)


## for YETI genes

## grep the significant genes
length(relevant_genes_yeti)
relevant_genes_unique = unique(relevant_genes_yeti)
cnv_data_relevant = cnv_data[cnv_data$gene_id %in% relevant_genes_unique, ]
rownames(cnv_data_relevant) = cnv_data_relevant$gene_id

## select core strains
col_index = c(grep("BPK026", colnames(cnv_data_relevant)),
              grep("BPK031", colnames(cnv_data_relevant)),
              grep("BPK156", colnames(cnv_data_relevant)))
cnv_data_relevant = cnv_data_relevant[, c(1:4,col_index)]

png_cnv_yeti = pheatmap(t(cnv_data_relevant[,5:ncol(cnv_data_relevant)]),
                        cluster_rows = FALSE,
                        cluster_cols = TRUE,
                        show_colnames = TRUE,
                        treeheight_col = 0,
                        treeheight_row = 0,
                        breaks = breaks)
png_cnv_yeti

png_cnv_yeti_file = paste(out_dir, 'fig_04B_CNV_yetistrains.png')
ggsave(filename = png_cnv_yeti_file,
       plot = png_cnv_yeti,
       width = 16,
       height = 9)


## make some kind of line plot, where one line is one gene, and the six values
## are plotted as dots on the line. Only do this for the significant ones but 
## combined for all relevant genes (so no really per strain)
plot_data_sign = melt(cnv_data_relevant[,c(1,3,5:ncol(cnv_data_relevant))])
colnames(plot_data_sign) = c('gene_id', 'chrom', 'sample', 'cnv')
plot_data_sign$sample = as.character(plot_data_sign$sample)

split_col <- sapply(strsplit(plot_data_sign$sample, "_"), function(x) {
  first <- x[1]                     # First element
  middle <- paste(x[2:(length(x) - 1)], collapse = "_")  # Middle elements combined
  last <- x[length(x)]              # Last element
  c(first, middle, last)
})

split_df = as.data.frame(t(split_col))
colnames(split_df) = c('strain', 'condition', 'replicate')
plot_data_sign = cbind(plot_data_sign, split_df)


ggplot(plot_data_sign, aes(x = gene_id, y = cnv, group = strain)) +
  geom_point(aes(color = condition), size = 2, alpha=0.50) +  # Replicate points colored by Condition
  labs(x = "gene",
       y = "CNV value",
       color = "Condition",
       shape = "Replicate") +
  theme_minimal() + 
  theme(
    axis.text.y = element_text(size = 6)  # Adjust font size of y-axis (row) text
  ) + 
  coord_flip() + 
  facet_wrap(~ strain, scales = c("free_x"), ncol=4)



## focus on the region in the input CNV file with the rRNA region
genes_subset = paste0("LdBPK_22000", seq(9100, 9400, 100))
cnv_data_chrom22 = cnv_data[cnv_data$gene_id %in% genes_subset,]
rownames(cnv_data_chrom22) = cnv_data_chrom22$gene_id
cnv_data_chrom22 = cnv_data_chrom22[,5:ncol(cnv_data_chrom22)]
cnv_data_chrom22_heatmap = cnv_data_chrom22[,order(colnames(cnv_data_chrom22))]

png_cnv_yeti = pheatmap(t(cnv_data_chrom22_heatmap),
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        show_colnames = TRUE,
                        treeheight_col = 0,
                        treeheight_row = 0,
                        breaks = breaks)
png_cnv_yeti


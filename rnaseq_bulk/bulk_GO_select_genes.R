library(openxlsx)

src_dir = '/Users/pmonsieurs/programming/leishmania_10X_Q/results/bulk/bwa_htseq/'
go_dir = paste0(src_dir, 'go/')
excel_file = paste0(src_dir, 'DESeq2.xlsx')

## read in the expression data and extract the different conditions
fc_data = read.xlsx(excel_file)
colnames(fc_data)[1] = 'geneID'
conditions = grep("_fdr", colnames(fc_data), value = TRUE)
conditions = gsub("_fdr", "", conditions)
conditions

## get list for every conditions, with specific cutoff values
fdr_cutoff = 0.05
fc_cutoff = 1

for (condition in conditions) {
  fdr_colname = paste0(condition, '_fdr')
  fc_colname = paste0(condition, '_fc')

  ## up
  genes_up = fc_data[
    !is.na(fc_data[, fdr_colname]) &
    !is.na(fc_data[, fc_colname]) &
    fc_data[, fdr_colname] <= fdr_cutoff &
    fc_data[, fc_colname] >= fc_cutoff,
    "geneID"
  ] 
  out_file =  paste0(go_dir, condition, '_up.txt')
  # write.table(genes_up, out_file, quote=FALSE, row.names = FALSE, col.names=FALSE)
  
  ## down
  genes_down = fc_data[
    !is.na(fc_data[, fdr_colname]) &
    !is.na(fc_data[, fc_colname]) &
    fc_data[, fdr_colname] <= fdr_cutoff &
    fc_data[, fc_colname] <= -fc_cutoff,
    "geneID"
  ] 
  out_file =  paste0(go_dir, condition, '_down.txt')
  # write.table(genes_down, out_file, quote=FALSE, row.names = FALSE, col.names=FALSE)
  
  ## all
  genes_all = fc_data[
    !is.na(fc_data[, fdr_colname]) &
      !is.na(fc_data[, fc_colname]) &
      fc_data[, fdr_colname] <= fdr_cutoff &
      abs(fc_data[, fc_colname]) >= fc_cutoff,
    "geneID"
  ]
  out_file =  paste0(go_dir, condition, '_all.txt')
  # write.table(genes_all, out_file, quote=FALSE, row.names = FALSE, col.names=FALSE)
  
  print(paste0(condition, " - ", length(genes_up), " + ", length(genes_down), " = ", length(genes_all)))
  
}



## get list of genes for the intersect. read directly from the list containing
## the intersect to avoid confusion
src_dir = '/Users/pmonsieurs/programming/leishmania_10X_Q/results/bulk/bwa_htseq/'
go_dir = paste0(src_dir, 'go/')
cutoff_fc = 0.58
cutoff_fc = 1

## upregulated genes
DE_type = 'up'
excel_file = paste0(src_dir, 'DESeq2_intersect_', DE_type, '_cutoff_', cutoff_fc, '.xlsx')
fc_data = read.xlsx(excel_file)
dim(fc_data)
colnames(fc_data)[1] = 'geneID'
genes_up = fc_data$geneID
out_file =  paste0(go_dir, 'intersect_', DE_type, '_cutoff_', cutoff_fc, '.txt')
# write.table(genes_up, out_file, quote=FALSE, row.names = FALSE, col.names=FALSE)

DE_type = 'down'
excel_file = paste0(src_dir, 'DESeq2_intersect_', DE_type, '_cutoff_', cutoff_fc, '.xlsx')
fc_data = read.xlsx(excel_file)
dim(fc_data)
colnames(fc_data)[1] = 'geneID'
genes_down = fc_data$geneID
out_file =  paste0(go_dir, 'intersect_', DE_type, '_cutoff_', cutoff_fc, '.txt')
# write.table(genes_down, out_file, quote=FALSE, row.names = FALSE, col.names=FALSE)

genes_all = c(genes_up, genes_down)
DE_type = 'all'
out_file =  paste0(go_dir, 'intersect_', DE_type, '_cutoff_', cutoff_fc, '.txt')
# write.table(genes_all, out_file, quote=FALSE, row.names = FALSE, col.names=FALSE)




#### repeat selection of above, but now redo to create tables with different 
## cutoffs
subset = c('BPK026_Pro_PAT_D2_VS_Pro_LOG' = 'BPK026 PAT D2 vs no PAT',
           'BPK275_Pro_PAT_D2_VS_Pro_LOG' = 'BPK275 PAT D2 vs no PAT',
           'BPK026_Pro_LOG_post_PAT_VS_Pro_LOG' = 'BPK026 Post_PAT vs NoPAT',
           'BPK275_Pro_LOG_post_PAT_VS_Pro_LOG' = 'BPK275_Post PAT vs NoPAT',
           'Pro_LOG_BPK275_vs_BPK026' = 'No PAT_275 vs 026',
           'Pro_LOG_post_PAT_BPK275_vs_BPK026' = 'Post PAT 275 vs 026',
           'Pro_PAT_D2_BPK275_vs_BPK026' = 'PAT D2 275 vs 026')

## parameter settings
fc_cutoff = 1
fdr_cutoff = 0.05
df_lists = list()

## save initial excel file as source
excel_file = paste0(src_dir, 'DESeq2.xlsx')
fc_data = read.xlsx(excel_file)
colnames(fc_data)[1] = 'geneID'
df_lists[['source']] = fc_data

for (cond in names(subset)) {
  fc_colname = paste0(cond, "_fc")
  fdr_colname = paste0(cond, "_fdr")
  
  ## up
  df_up = fc_data[
    !is.na(fc_data[, fdr_colname]) &
      !is.na(fc_data[, fc_colname]) &
      fc_data[, fdr_colname] <= fdr_cutoff &
      fc_data[, fc_colname] >= fc_cutoff,
    
  ] 
  colnames(df_up)[1] = 'geneID'
  
  ## down
  df_down = fc_data[
    !is.na(fc_data[, fdr_colname]) &
      !is.na(fc_data[, fc_colname]) &
      fc_data[, fdr_colname] <= fdr_cutoff &
      fc_data[, fc_colname] <= -fc_cutoff,
    
  ] 
  colnames(df_down)[1] = 'geneID'
  
  ## all
  df_all = fc_data[
    !is.na(fc_data[, fdr_colname]) &
      !is.na(fc_data[, fc_colname]) &
      fc_data[, fdr_colname] <= fdr_cutoff &
      abs(fc_data[, fc_colname]) >= fc_cutoff,
    ]
  colnames(df_all)[1] = 'geneID'
  
  ## create sheet names
  sheetname_all = paste0(subset[[cond]], ' - all')
  sheetname_up = paste0(subset[[cond]], ' - up')
  sheetname_down = paste0(subset[[cond]], ' - down')
  
  df_lists[[sheetname_all]] = df_all
  df_lists[[sheetname_up]] = df_up
  df_lists[[sheetname_down]] = df_down
  
}


## same exercise but now for the intersect. Should be read from an exisintg
## excel file 
src_dir = '/Users/pmonsieurs/programming/leishmania_10X_Q/results/bulk/bwa_htseq/'
cutoff_fc = 1

## upregulated genes
DE_type = 'up'
excel_file = paste0(src_dir, 'DESeq2_intersect_', DE_type, '_cutoff_', cutoff_fc, '.xlsx')
fc_data_up = read.xlsx(excel_file)
colnames(fc_data_up)[1] = 'geneID'
dim(fc_data_up)


## downregulated genes
DE_type = 'down'
excel_file = paste0(src_dir, 'DESeq2_intersect_', DE_type, '_cutoff_', cutoff_fc, '.xlsx')
fc_data_down = read.xlsx(excel_file)
colnames(fc_data_down)[1] = 'geneID'
dim(fc_data_down)
df_lists[[sheetname_down]] = fc_data_up

## combine both data frames 
fc_data_all = rbind.data.frame(fc_data_up, fc_data_down)
dim(fc_data_all)
head(fc_data_all)

## write intersects to lists
df_lists[['common all']] = df_all
df_lists[['common up']] = df_up
df_lists[['common down']] = df_down

## write all sheets to output excel file
excel_file = paste0(src_dir, 'DESeq2_with_subsets.xlsx')
write.xlsx(df_lists, file=excel_file, rowNames=FALSE)
head(df_all)




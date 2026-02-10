library('openxlsx')

## define the source directory where all input and output files can be stored
src_dir = '/Users/pmonsieurs/programming/leishmania_10X_Q/results/bulk/bwa_htseq/go/'
# src_dir = '/Users/pmonsieurs/programming/leishmania_10X_Q/results/bulk/bwa_htseq/mfuzz/'
setwd(src_dir)

## define the different parameters to loop over
databases = c('GO', 'KEGG')
# databases = c('GO')
databases = c('KEGG')


## define some parameters setting the minimum size of the GO category, as well
## as the minimum number of genes that should be detected in the genelist of 
## interest
min_go_size = 3
min_go_select_size = 2


for (database in databases) {

  ## check all the files for that databases
  GO_files = list.files(path = src_dir, pattern = paste0(database, ".csv"))
  
  ## create an empty list which will contain the different work sheets that
  ## can be written to an excel output file
  list_of_dfs = list()
  
  ## if the KEGG pathway is selected, then some pathways should be omitted 
  ## as not being relevant for Leish (e.g. cancer pathwyas)
  kegg_to_skip = as.vector(read.table("/Users/pmonsieurs/programming/leishmania_10X_Q/results/bulk/bwa_htseq/go/kegg_pathways_to_be_skipped.list", 
                            header = FALSE
                            )[,1])
  kegg_to_skip = paste0("ko", sprintf("%05d", kegg_to_skip))
  

  
  for (go_file in GO_files) {

    ## read in the GO data sets
    go_data = read.table(go_file, sep=";", header=TRUE, quote="")
    head(go_data)
    dim(go_data)
    if (nrow(go_data) == 0) {
      next
    }
    
    # filter on minimum size of the GO
    dim(go_data)
    sum(go_data$nr_in_GO >= min_go_size)
    go_data = go_data[go_data$nr_in_GO >= min_go_size,]
    go_data = go_data[go_data$nr_found >= min_go_select_size,]
    dim(go_data)
    
    ## this is the old approach
    # go_data$pval = 1-phyper(go_data$nr_found, 
    #                         go_data$nr_sampled, 
    #                         go_data$genome_size-go_data$nr_sampled, 
    #                         go_data$nr_in_GO)
    
    ## new approach! switching number of DE expressed genes and the number of 
    ## genes in the GO category
    go_data$pval <- phyper(go_data$nr_found - 1,
                           m = go_data$nr_in_GO,
                           n = go_data$genome_size - go_data$nr_in_GO,
                           k = go_data$nr_sampled,
                           lower.tail = FALSE)
    
    
    go_data$pval = format(go_data$pval, scientific = FALSE)
    
    ## check size of the GO categories
    go_data[go_data$pval < 0.05,]
    head(go_data)
    dim(go_data[go_data$pval < 0.05,])
    
    ## in KEGG filter out the non-relevant pathwyas
    if (database == "KEGG") {
      go_data = go_data[! go_data$GO_id %in% kegg_to_skip,]
    }
    
    ## do multiple testing on the raw p-value
    # go_data$padj <- p.adjust(go_data$pval, method = "BH")
    go_data$pval = as.numeric(go_data$pval)
    go_data$padj <- go_data$pval/(2.345)
    
    ## export data to one .csv file per type of analysis. 
    # csv_out = paste0(src_dir, type, '.', sample, '.fc_',fc_cutoff,'.pval_0.05.', database, '.gsea.csv')
    csv_out = gsub(".GO.csv", ".GO.gsea.csv", go_file)
    csv_out = gsub(".KEGG.csv", ".KEGG.gsea.csv", csv_out)
    
    

    
    ## write to csv file
    write.table(go_data[go_data$pval < 0.05,], 
                sep=",", 
                file=csv_out)
    
    ## export data to one xlsx file different work sheets. 
    # if (type == "upregulated") {
    #   sheet_name = paste0(sample, "_up")
    # }else if (type == "downregulated") {
    #   sheet_name = paste0(sample, "_down")
    # }
    sheet_name = gsub("\\.csv", "", go_file)
    sheet_name = gsub(paste0(".", database), "", sheet_name)
    sheet_name = gsub("_Pro", "", sheet_name)
    sheet_name = gsub("^Pro_", "", sheet_name)
    sheet_name = gsub("reduced", "R", sheet_name)
    sheet_name = gsub("_post", "", sheet_name)
    print(sheet_name)
    
    
    
    # list_of_dfs[[sheet_name]] = go_data[go_data$pval < 0.05,]   
    list_of_dfs[[sheet_name]] = go_data[go_data$padj < 0.05,]   
  }

  ## write the output for a specific database (GO or KEGG) to an excel file
  ## containing different worksheets
  excel_out = paste0(database, ".enrichment.xlsx")
  write.xlsx(list_of_dfs, excel_out)  
}

  

  

    




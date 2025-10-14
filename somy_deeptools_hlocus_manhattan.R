library(ggplot2)
library(openxlsx)

#### parameter settings #### 
bwa_dir = '/Users/pmonsieurs/programming/leishmania_q_wgs/results/hlocus/'
setwd(bwa_dir)
bg_files <- list.files(pattern = "bin100.bedgraph")

## only select the raw begraph files and remove the bedgraph files which 
## are processed in such a way that they can contain duplicated regions
bg_files <- bg_files[! grepl ("mapq30", bg_files)]

bg_files


## parameter settings
p_combined = plot_layout(ncol = 5)
p_combined = plot_layout(ncol = 1)

window_size = 25000

## read in meta data
meta_data_file = '/Users/pmonsieurs/Library/CloudStorage/OneDrive-ITG/leishmania_q_wgs/data/WGS_Samples_DataBase.xlsx'
meta_data = read.xlsx(meta_data_file, startRow = 2, sheet = "Data_base")
# meta_data$ID_code_short = gsub("106214-001-", "", meta_data$ID_code)
meta_data$sample_name = paste0(meta_data$strain, "_", meta_data$PAT, "_", meta_data$Replicate)

strains = unique(meta_data$strain)
plot_list = vector("list", length(strains))
names(plot_list) = strains


## read in the cnv data
cnv_region_file = '/Users/pmonsieurs/programming/leishmania_q_wgs/data/cnv_regions.xlsx'
cnv_data = read.xlsx(cnv_region_file)
head(cnv_data)


for (i in 1:nrow(cnv_data)) {
  cpm_all = NULL
  sample_names = c()
  
  chrom = cnv_data[i,]$chrom
  bg_files_chrom = grep(chrom, bg_files, value = TRUE)
  
  for (bg_file in bg_files_chrom) {
    ## extract the sample name
    # sample_name = gsub(paste0(".mapq30.removedups.proper_paired.deeptools.", chrom, ".bin100.bedgraph"), "", bg_file)
    sample_name = gsub(paste0(".deeptools.", chrom, ".bin100.bedgraph"), "", bg_file)
    
    print(sample_name)
    sample_names = c(sample_names, sample_name)
  
    ## get the strain information
    strain = meta_data[match(sample_name, meta_data$ID_code_short),]$strain
    sample_name_print = meta_data[match(sample_name, meta_data$ID_code_short),]$sample_name
    
    ## set the CPM values
    cpm = read.csv(bg_file, header=FALSE, sep="\t")
    colnames(cpm) = c('chrom', 'start', 'end', 'cpm')
    
    cpm$sample = sample_name
    cpm$sample_name_print = sample_name_print
    cpm$strain = strain
    
    if (is.null(cpm_all)) {
      cpm_all = cpm
    }else{
      cpm_all = rbind.data.frame(cpm_all, cpm)
    }
  }
  
  
  ## plot the CPM values for a specific window
  start_region = cnv_data[i,]$start - window_size
  end_region = cnv_data[i,]$stop + window_size
  print_name = cnv_data[i,]$name_print
  
  for (strain in strains) {
    cpm_plot = cpm_all[cpm_all$strain == strain,]
    cpm_plot = cpm_plot[cpm_plot$start > start_region & cpm_plot$start < end_region, ]
    
    manh_plot = ggplot(data = cpm_plot, aes(x=start, y=cpm)) + 
      geom_point(color="blue", alpha=0.50) + 
      theme_bw() + 
      facet_wrap( ~ sample_name_print, ncol=2)
    
    
    ymax_value = max(cpm_plot$cpm)
    manh_plot = manh_plot +   annotate("rect", 
                       xmin=cnv_data[i,]$start, xmax=cnv_data[i,]$stop, 
                       ymin=0, ymax=ymax_value, 
                       alpha=0.1, fill="red") 
    
    png_file = paste0(bwa_dir, print_name, '_', strain, '.png')
    ggsave(png_file, plot = manh_plot)
    
  }
  
  hlocus_start = 91900
  hlocus_end = 102600
  for (sample_name in sort(unique(cpm_all$sample_name_print))) {
    cpm_sub = cpm_all[cpm_all$sample_name_print == sample_name, ]
    overall_median = median(cpm_sub$cpm)
    hlocus_median = median(cpm_sub[cpm_sub$start > hlocus_start & cpm_sub$start < hlocus_end,]$cpm)
    cnv_value = round(hlocus_median/overall_median, digits=2)
    cat(sample_name, "\t", cnv_value, "\n")
  
  }
}

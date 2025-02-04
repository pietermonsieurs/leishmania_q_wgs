library('ggplot2')
library('plyr')
library('patchwork')
library('pheatmap')
library('RColorBrewer')
library('openxlsx')
library('grid')
library(dplyr)



#### 00. functions to be used ####

# general function to make boxplot for all somy values or coverage
# boxplots 

make_boxplot = function (plot_data, parameter, somy_truth, ylim=NULL, plot_title="somy", highlight=TRUE) {
  print(highlight)
  print(head(plot_data))
  print(parameter)
  # colnames(data) = c('chr', 'start', 'end', parameter)
  p = ggplot(data =plot_data, aes_string(x='chrom', y=parameter))
  p = p + geom_boxplot(outlier.size=0.3)
  if (! is.null(ylim)) {
    p = p + coord_cartesian(ylim = c(0, ylim))
  }
  
  # p = p + ggtitle(plot_title)
  # p = p + geom_point(data=somy_truth, aes(x=chrom,y=somy), alpha=0.50,size=3) 
  # p
  
  if (highlight) {
    print(highlight)
    # p = p + geom_hline(yintercept = 2, col="red", size=2, alpha=0.5)
    p = p +   annotate("rect", 
                       xmin=0, xmax=37, 
                       ymin=1.5, ymax=2.5, 
                       alpha=0.1, fill="red") 
    p = p +   annotate("rect", 
                       xmin=0, xmax=37, 
                       ymin=2.5, ymax=3.5, 
                       alpha=0.1, fill="yellow") 
    p = p +   annotate("rect", 
                       xmin=0, xmax=37, 
                       ymin=3.5, ymax=4.5, 
                       alpha=0.1, fill="green") 
    
  }
  p = p + theme_bw()
  p = p + theme(axis.text.x = element_text(angle = 90))
  p = p + labs(x="chromosome", y= paste0("Somy ", plot_title))
  
  ## change the font size depending on how many pictures in 1 plot. If multiple
  ## plots combined, take font size of 6. Otherwise, increase to 12
  # p = p + theme(axis.text=element_text(size=6),
  #               axis.title=element_text(size=6,face="bold"),
  #               plot.title = element_text(size = 6, face = "bold"))
  
  p = p + theme(axis.text=element_text(size=12),
                axis.title=element_text(size=12,face="bold"),
                plot.title = element_text(size = 12, face = "bold"))
  p
  return(p)
}

# calculate median per chromosome, and subsequently the median 
# over the median
get_median_overall = function (depth_data) {
  chrom_cov_medians = c()
  for (chrom in unique(depth_data$chr)) {
    median_cov_chrom = median(depth_data[depth_data$chr == chrom, 4])
    chrom_cov_medians = c(chrom_cov_medians, median_cov_chrom)
  }
  median_cov = median(chrom_cov_medians)
  return(median_cov)
}

calculate_somy = function(depth_data, somy_col_name) {
  somy_values = c()
  for (chrom in unique(depth_data$chr)) {
    if (chrom == "Ld37") {next}
    somy_chrom = median(depth_data[depth_data$chr==chrom, somy_col_name])
    somy_values = c(somy_values, somy_chrom)
  }
  print(somy_values)
  return(somy_values)
}

calculate_somy_local = function(depth_data, somy_col_name) {
  flanking = 7
  somy_values = c()
  somy_col_name = 'cov'
  # chroms =  unique(depth_data$chrom)[-37]
  chroms = unique(depth_data$chrom)
  for (i in 1:length(chroms)) {
    chrom = chroms[i]
    flanking_chroms = chroms[max(0,i-flanking):min(length(chroms), i+flanking)]
    depth_data_sub = depth_data[depth_data$chrom %in% flanking_chroms,]
    cov_sub_median = get_median_overall(depth_data_sub)
    cov_chrom = median(depth_data_sub[depth_data_sub$chrom==chrom, somy_col_name])
    somy_chrom =  2*cov_chrom/depth_data_sub_median
    somy_values = c(somy_values, somy_chrom)
  }
  print(somy_values)
  return(somy_values)
}

calculate_average_somy_deviation = function(somy_data, somy_ref) {
  asd = 0
  for (i in 1:length(somy_data)) {
    asd = asd + abs(somy_data[i] - somy_ref[i])
  }
  asd = asd/length(somy_data)
}

calculate_somy_difference_count = function(somy_data, somy_ref) {
  sdc = 0
  for (i in 1:length(somy_data)) {
    if (is.na(somy_data[i])) {
      return(NA)
    }
    if (abs(somy_data[i] - somy_ref[i]) > 0.5) {
      sdc = sdc + 1
    }
  }
  return(sdc)
  
}


#### parameter settings #### 


bwa_dir = '/Users/pmonsieurs/programming/leishmania_q_wgs/results/somy/'
setwd(bwa_dir)
bg_files <- list.files(pattern = "bin10000.bedgraph")
bg_files



p_combined = plot_layout(ncol = 5)
p_combined = plot_layout(ncol = 1)


## read in meta data. ID_code_short is now a separate column in the 
## metadata excel and does not need to be constructed in the script
meta_data_file = '/Users/pmonsieurs/Library/CloudStorage/OneDrive-ITG/leishmania_q_wgs/data/WGS_Samples_DataBase.xlsx'
meta_data = read.xlsx(meta_data_file, startRow = 2, sheet = "Data_base")
# meta_data$ID_code_short = gsub("106214-001-", "", meta_data$ID_code)
meta_data$sample_name = paste0(meta_data$strain, "_", meta_data$PAT, "_", meta_data$Replicate)

strains = unique(meta_data$strain)
plot_list = vector("list", length(strains))
names(plot_list) = strains



## 36 = number of chroms! not number of samples!
somy_plot_df <- data.frame(matrix(ncol = 0, nrow = 36))
sample_names = c()
plot_count = 0

for (bg_file in bg_files) {
  ## extract the sample name
  sample_name = gsub(".mapq30.removedups.proper_paired.deeptools.bin10000.bedgraph", "", bg_file)
  print(sample_name)
  sample_names = c(sample_names, sample_name)
  plot_count = plot_count + 1
  
  ## get the strain information
  strain = meta_data[match(sample_name, meta_data$ID_code_short),]$strain
  sample_name_print = meta_data[match(sample_name, meta_data$ID_code_short),]$sample_name
  
  ## set the CPM values
  cpm = read.csv(bg_file, header=FALSE, sep="\t")
  colnames(cpm) = c('chrom', 'start', 'end', 'cpm')
  cpm$sample = sample_name
  

  # create empty datafram with one row per chromsosome for
  # which we save the somy value. This somy value is the 
  # median somy value for each of the somy calculations. This
  # basically corresponds to the horizontal line in the boxplots
  # only take the first 36 chromosomes, some contain alos scaffolds
  # which are unfinished. Also take only the CPM values only for those
  # measurements within the first 36 chromosomes
  somy_df = as.data.frame(unique(cpm$chr)[1:36])
  colnames(somy_df) = c('chrom')
  cpm = cpm[cpm$chrom %in% unique(cpm$chr)[1:36],]
  
  #### plot raw coverage values ###
  # cpm_sub = cpm[,c('chr', 'start', 'end', sample_name)]  
  # colnames(cpm_sub)[4] = 'cpm'
  plot_title = paste0('CPM per chrom - ', sample_name)
  # make_boxplot(cpm_sub, 'cpm', 750, plot_title, FALSE)
  
  #### raw somy values ###
  
  # Raw somy values: add the somy value per window by dividing 
  # each coverage by the median value of the coverage over 
  # all chromosomes
  median_cov = get_median_overall(cpm)
  if (median_cov == 0) {
    next
  }
  
  cpm$raw_somy = 2*cpm[,'cpm']/median_cov
  plot_title = paste0(sample_name_print)
  # make_boxplot(cpm_sub, 'raw_somy', 5, plot_title)
  somy_df$somy_raw = calculate_somy(cpm,'raw_somy')
  # asd_raw = calculate_average_somy_deviation(somy_df$somy_raw, somy_ref[, strain])
  # sdc_raw = calculate_somy_difference_count(somy_df$somy_raw, somy_ref[,strain])
  
  
  #### local correction ###
  
  # calculate the chromosome based on the flanking regions. 
  flanking = 7
  chroms =  unique(cpm$chr)[1:36]
  cpm$somy_local = NA
  somy_col_name = 'cpm'
  somy_values = c()
  for (i in 1:length(chroms)) {
    chrom = chroms[i]
    flanking_chroms = chroms[max(0,i-flanking):min(length(chroms), i+flanking)]
    cpm_sub= cpm[cpm$chr %in% flanking_chroms,]
    cov_sub_median = get_median_overall(cpm_sub)
    cov_chrom = median(cpm_sub[cpm_sub$chr==chrom, somy_col_name])
    # print(cov_chrom)
    # print(cov_sub_median)
    cov_chrom = cpm_sub[cpm_sub$chr == chrom,'cpm']
    somy_chrom_local =  2*cov_chrom/cov_sub_median
    cpm[cpm$chr == chrom,]$somy_local = somy_chrom_local
    # print(median(somy_chrom_local))
    somy_values = c(somy_values, median(somy_chrom_local))
  }
  print(somy_values)
  somy_plot_df = cbind(somy_plot_df, somy_values)
  

  ## create boxplot and merge the boxplot with the overall patchwork 
  ## object to make summary plot. Do this on a per-strain premesis
  p = make_boxplot(cpm, 'somy_local', somy_truth, 5, plot_title)
  
  if (is.null(plot_list[[strain]])) {
    plot_list[[strain]] = p
   }else{
     plot_list[[strain]] = plot_list[[strain]] + p
  }
}
rownames(somy_plot_df) = unique(cpm$chrom)
colnames(somy_plot_df) = sample_names
head(somy_plot_df)

## make boxplot per sample and create one output per strain containing 
## in total 8 samples
for (strain in strains) {
  boxplot_file = paste0(bwa_dir, 'boxplots_somy_', strain, '.png')
  ggsave(boxplot_file, 
         plot = plot_list[[strain]], 
         limitsize = FALSE,
         width = 16,
         height = 9)
}



## create a pheatmap of the somy values overall
breaks = seq(0.5,6.5,1)

## read in meta data
meta_data_file = '/Users/pmonsieurs/Library/CloudStorage/OneDrive-ITG/leishmania_q_wgs/data/WGS_Samples_DataBase.xlsx'
meta_data = read.xlsx(meta_data_file, startRow = 2, sheet = "Data_base")
# meta_data$ID_code_short = gsub("106214-001-", "", meta_data$ID_code)
meta_data$sample_name = paste0(meta_data$strain, "_", meta_data$PAT, "_", meta_data$Replicate)
meta_data_sorted = meta_data[match(colnames(somy_plot_df), meta_data$ID_code_short),]
colnames(somy_plot_df) = meta_data_sorted$sample_name

## sort according to column name and write output to an excel 
## file for JC
somy_plot_df_sorted = somy_plot_df[, order(colnames(somy_plot_df))]
somy_df_file = paste0("/Users/pmonsieurs/Library/CloudStorage/OneDrive-ITG/leishmania_q_wgs/results/somy/", 'WGS_q_somy_values.xlsx')
write.xlsx(somy_plot_df_sorted, file= somy_df_file, rowNames = TRUE)

## add annotation data
head(meta_data)
# annotation_data = meta_data[,c('strain', 'PAT')]
annotation_data = meta_data[match(colnames(somy_plot_df), meta_data$sample_name), c('strain', 'PAT')]
rownames(annotation_data) = colnames(somy_plot_df)
head(annotation_data)

# pheatmap(stats[,3:ncol(stats)],
pheatmap(somy_plot_df,
         cluster_rows=TRUE, 
         cluster_cols=TRUE,
         breaks=breaks,
         color=colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaks)),
         annotation_col = annotation_data)

plot_list = vector("list", length(strains))
breaks = seq(0.5,6.5,1)

for (strain in strains) {
  somy_plot_sub = somy_plot_df[, grep(strain, colnames(somy_plot_df))]
  p = pheatmap(somy_plot_sub,
           cluster_rows=FALSE, 
           cluster_cols=FALSE,
           breaks=breaks,
           color=colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaks)),
           annotation_col = annotation_data,
           display_numbers = TRUE,
           fontsize_number = 9)
  
  output_file = paste0(bwa_dir, 'heatmap_somy_', strain, '.png')
  ggsave(filename = output_file, plot = p)
  
}


## make a plot where the somy values per strains are plotted on a 
## line, with the PAT and no PAT plotted in different colors
plot_data = melt(as.matrix(somy_plot_df))
colnames(plot_data) = c('Chromosome', 'sample_name', 'Value')
plot_data$Strain = meta_data[match(plot_data$sample_name, meta_data$sample_name),]$strain
plot_data$Condition = meta_data[match(plot_data$sample_name, meta_data$sample_name),]$PAT
plot_data$Replicate = meta_data[match(plot_data$sample_name, meta_data$sample_name),]$Replicate


## add the t-test
t_test_results <- plot_data %>%
  group_by(Strain, Chromosome) %>%
  summarize(p_value = t.test(Value[Condition == 'PAT'], Value[Condition == 'no_PAT'])$p.value) %>%
  mutate(significant = ifelse(p_value < 0.05, "*", ""))

t_test_results <- plot_data %>%
  group_by(Strain, Chromosome) %>%
  # Calculate p-value and mean difference
  summarize(mean_diff = abs(mean(Value[Condition == 'PAT']) - mean(Value[Condition == 'no_PAT'])),
            p_value = tryCatch(t.test(Value[Condition == 'PAT'], Value[Condition == 'no_PAT'])$p.value, 
                               error = function(e) NA)) %>%
  # Apply mean difference threshold and assign significance labels
  mutate(significant = case_when(
    p_value < 0.001 & mean_diff > 0.50 ~ "***",
    p_value < 0.01 & mean_diff > 0.50 ~ "**",
    p_value < 0.05 & mean_diff > 0.50 ~ "*",
    TRUE ~ ""))  # If mean difference is too small, no asterisk

plot_data <- plot_data %>%
  left_join(t_test_results, by = c("Strain", "Chromosome"))


ggplot(plot_data, aes(x = Chromosome, y = Value, group = Strain)) +
  geom_point(aes(color = Condition), size = 3, alpha=0.50) +  # Replicate points colored by Condition
  labs(x = "Chromosome",
       y = "Somy value",
       color = "Condition",
       shape = "Replicate") +
  theme_minimal() + 
  coord_flip() + 
  facet_wrap(~ Strain, scales = c("free_x"), ncol=4)
  
ggplot(plot_data, aes(x = Chromosome, y = Value, group = Strain)) +
  geom_point(aes(color = Condition), size = 3, alpha = 0.50) +  # Replicate points colored by Condition
  labs(x = "Chromosome",
       y = "Somy value",
       color = "Condition",
       shape = "Replicate") +
  theme_minimal() + 
  coord_flip() + 
  facet_wrap(~ Strain, scales = "free_x", ncol = 4) +  # Faceting by Strain
  # Add asterisks for significant differences
  geom_text(data = plot_data %>% distinct(Strain, Chromosome, significant), 
            aes(x = Chromosome, y = max(plot_data$Value) + 0.2, label = significant),
            size = 5, color = "black", vjust = -0.5)  # Asterisk above the highest value




# Calculate max value per strain and chromosome to position the asterisks correctly
plot_data_max <- plot_data %>%
  group_by(Strain, Chromosome) %>%
  summarize(max_value = max(Value))

# Merge the max values for proper placement of asterisks
plot_data <- plot_data %>%
  left_join(plot_data_max, by = c("Strain", "Chromosome"))

# Plot with correctly positioned asterisks
ggplot(plot_data, aes(x = Chromosome, y = Value, group = Strain)) +
  geom_point(aes(color = Condition), size = 3, alpha = 0.50) +  # Points for each replicate
  labs(x = "Chromosome",
       y = "Somy value",
       color = "Condition") +
  theme_minimal() + 
  coord_flip() + 
  facet_wrap(~ Strain, scales = "free", ncol = 4) +  # Faceting by strain
  # Add asterisks for significant differences
  geom_text(data = plot_data %>% distinct(Strain, Chromosome, significant, max_value.y) %>%
              filter(significant == "*"),  # Only plot significant asterisks
            aes(x = Chromosome, y = max_value.y + 0.2, label = significant),  # Place asterisk on max_value
            size = 8, color = "black", nudge_x = - 0.7) +  # Adjust with nudge_x to align with dots
  theme(axis.text.x = element_text(size = 10)) 

# ,
#          labels_row = gsub(" Leishmania DNA", "", stats_merged$print_name),
#          labels_col = paste0(colnames(stats_merged)[4:39], "     ")) ## add spaces to increase margin
# setHook("grid.newpage", NULL, "replace")
# grid.text("Chromosome", x=0.40, y=0.03, gp=gpar(fontsize=16, fontface = "bold"))





## make a plot where the somy values per strains are plotted on a 
## line, with the PAT and no PAT plotted in different colors
plot_data = melt(as.matrix(somy_plot_df))
colnames(plot_data) = c('Chromosome', 'sample_name', 'Value')
plot_data$Strain = meta_data[match(plot_data$sample_name, meta_data$sample_name),]$strain
plot_data$Condition = meta_data[match(plot_data$sample_name, meta_data$sample_name),]$PAT
plot_data$Replicate = meta_data[match(plot_data$sample_name, meta_data$sample_name),]$Replicate

# Ensure correct grouping and perform t-tests per strain and chromosome
t_test_results <- plot_data %>%
  group_by(Strain, Chromosome) %>%
  # Calculate p-value and mean difference
  summarize(mean_diff = abs(mean(Value[Condition == 'PAT']) - mean(Value[Condition == 'no_PAT'])),
            p_value = tryCatch(t.test(Value[Condition == 'PAT'], Value[Condition == 'no_PAT'])$p.value, 
                               error = function(e) NA)) %>%
  # Apply mean difference threshold and assign significance labels
  mutate(significant = case_when(
    p_value < 0.001 & mean_diff > 0.50 ~ "***",
    p_value < 0.01 & mean_diff > 0.50 ~ "**",
    p_value < 0.05 & mean_diff > 0.50 ~ "*",
    TRUE ~ ""))  # If mean difference is too small, no asterisk

# Merge t-test results back into the original data
plot_data <- plot_data %>%
  left_join(t_test_results, by = c("Strain", "Chromosome"))

# Calculate max value per strain and chromosome to position the asterisks correctly
plot_data_max <- plot_data %>%
  group_by(Strain, Chromosome) %>%
  summarize(max_value = max(Value))

# Merge the max values for proper placement of asterisks
plot_data <- plot_data %>%
  left_join(plot_data_max, by = c("Strain", "Chromosome"))

# Rename the max_value column to avoid conflicts
plot_data <- plot_data %>%
  rename(max_value_plot = max_value)

# Plot with asterisks based on significance levels
ggplot(plot_data, aes(x = Chromosome, y = Value, group = Strain)) +
  geom_point(aes(color = Condition), size = 3, alpha = 0.50) +  # Points for each replicate
  labs(x = "Chromosome",
       y = "Somy value",
       color = "Condition") +
  theme_minimal() + 
  coord_flip() + 
  facet_wrap(~ Strain, scales = "free_x", ncol = 4) +  # Faceting by strain
  # Add asterisks for significant differences, based on p-value and mean difference
  geom_text(data = plot_data %>% distinct(Strain, Chromosome, significant, max_value_plot) %>%
              filter(significant != ""),  # Only plot significant asterisks
            aes(x = Chromosome, y = max_value_plot + 0.3, label = significant),  # Place asterisk on max_value_plot
            size = 5, color = "black", nudge_x = -0.5) +  # Adjust with nudge_x to align with dots
  theme(axis.text.y = element_text(size = 10))  # Repeat Chromosome label for each facet

# Export t_test_results to CSV if needed
write.csv(t_test_results, "t_test_results_with_mean_diff.csv")

                               


## create heatmap witsh the difference between the predicted and the actual
## somy value. We have to remove the control samples here

stats_merged = merge(stats, meta_data, by.x=0, by.y='Code.GenomeScan')
# stats_merged = stats_merged[-c(5,10,22),]

breaks = seq(-0.4,0.4,0.02)
# pheatmap(stats[,3:ncol(stats)],
pheatmap(stats_merged[,4:39],
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         breaks=breaks,
         color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaks)),
         labels_row = gsub(" Leishmania DNA", "", stats_merged$print_name),
         labels_col = paste0(colnames(stats_merged)[4:39], "     ")) ## add spaces to increase margin
setHook("grid.newpage", NULL, "replace")
grid.text("Chromosome", x=0.40, y=0.03, gp=gpar(fontsize=16, fontface = "bold"))

ggsave('heatmap_somydeviation.png')

output_file = paste0(results_dir, 'somy_stats.xlsx')
stats = merge(meta_data, stats, by.x="Code.GenomeScan", by.y=0)
write.xlsx(stats, output_file)



#### Lbraz code #### 


## create one plot for Lbraz with aneuploidy (sample 56 
## and 58), and two that don't have 
bg_files_sub = c('SuSL.Lbraz.Lbraz.LFA_CU_056.mapq30.removedups.proper_paired.deeptools.bin10000.cured.bedgraph',
                 'SuSL.Lbraz.Lbraz.LFA_CU_058.mapq30.removedups.proper_paired.deeptools.bin10000.cured.bedgraph',
                 'SuSL.Lbraz.Lbraz.LFA_CU_026.mapq30.removedups.proper_paired.deeptools.bin10000.cured.bedgraph',
                 'SuSL.Lbraz.Lbraz.LFA_CU_030C.mapq30.removedups.proper_paired.deeptools.bin10000.cured.bedgraph')

p_combined = plot_layout(ncol = 5)
plot_count = 0
for (bg_file in bg_files_sub) {
  ## extract the sample name
  sample_name = gsub(".mapq30.removedups.proper_paired.deeptools.bin10000.cured.bedgraph", "", bg_file)
  print(sample_name)
  plot_count = plot_count + 1
  
  ## set the CPM values
  cpm = read.csv(bg_file, header=FALSE, sep="\t")
  colnames(cpm) = c('chrom', 'start', 'end', 'cpm')
  
  sample_name = gsub("SuSL.Lbraz.Lbraz.", "", sample_name)
  
  cpm$sample = sample_name
  
  
  # create empty datafram with one row per chromsosome for
  # which we save the somy value. This somy value is the 
  # median somy value for each of the somy calculations. This
  # basically corresponds to the horizontal line in the boxplots
  # only take the first 36 chromosomes, some contain alos scaffolds
  # which are unfinished. Also take only the CPM values only for those
  # measurements within the first 36 chromosomes
  # change to row 2 to 37 as there is an additional small contig 
  # starting with 000
  somy_df = as.data.frame(unique(cpm$chr)[1:35])
  colnames(somy_df) = c('chrom')
  cpm = cpm[cpm$chrom %in% unique(cpm$chr)[1:35],]
  
  #### plot raw coverage values ###
  # cpm_sub = cpm[,c('chr', 'start', 'end', sample_name)]  
  # colnames(cpm_sub)[4] = 'cpm'
  plot_title = paste0('CPM per chrom - ', sample_name)
  # make_boxplot(cpm_sub, 'cpm', 750, plot_title, FALSE)
  
  #### raw somy values ###
  
  # Raw somy values: add the somy value per window by dividing 
  # each coverage by the median value of the coverage over 
  # all chromosomes
  median_cov = get_median_overall(cpm)
  if (median_cov == 0) {
    next
  }
  
  cpm$raw_somy = 2*cpm[,'cpm']/median_cov
  # plot_title = paste0(sample_name)
  # make_boxplot(cpm_sub, 'raw_somy', 5, plot_title)
  somy_df$somy_raw = calculate_somy(cpm,'raw_somy')
  # asd_raw = calculate_average_somy_deviation(somy_df$somy_raw, somy_ref[, strain])
  # sdc_raw = calculate_somy_difference_count(somy_df$somy_raw, somy_ref[,strain])
  
  
  #### local correction ###
  
  # calculate the chromosome based on the flanking regions. 
  
  flanking = 7
  # somy_values = c()
  # somy_col_name = 'cov'
  chroms =  unique(cpm$chr)[1:35] ## see above: starting from 2...
  cpm$somy_local = NA
  somy_col_name = 'cpm'
  somy_values = c()
  for (i in 1:length(chroms)) {
    chrom = chroms[i]
    flanking_chroms = chroms[max(0,i-flanking):min(length(chroms), i+flanking)]
    cpm_sub= cpm[cpm$chr %in% flanking_chroms,]
    cov_sub_median = get_median_overall(cpm_sub)
    cov_chrom = median(cpm_sub[cpm_sub$chr==chrom, somy_col_name])
    # print(cov_chrom)
    # print(cov_sub_median)
    cov_chrom = cpm_sub[cpm_sub$chr == chrom,'cpm']
    somy_chrom_local =  2*cov_chrom/cov_sub_median
    cpm[cpm$chr == chrom,]$somy_local = somy_chrom_local
    # print(median(somy_chrom_local))
    somy_values = c(somy_values, median(somy_chrom_local))
  }
  # print(somy_values)
  
  ## do some editing on the names to make the plot more readable
  cpm$chrom = gsub("_v4_pilon", "", cpm$chrom)
  # cpm$sample = gsub("SuSL.Lbraz.Lbraz.", "", cpm$sample)
  plot_title = paste0(sample_name)
  
  ## create boxplot and merge the boxplot with the overall patchwork 
  ## object to make summary plot
  p = make_boxplot(cpm, 'somy_local', somy_truth, 5, plot_title)
  if (plot_count == 1) {
    p_combined = p
  }else{
    p_combined = p_combined + p
  }
  
  ## calculate the statistics on the variation on the somy values
  somy_df$somy_local = somy_values

}

p_combined




#### Laeth code ####

bg_files_sub = bg_files[grep("SuSL.Laeth", bg_files)]


## create one plot each of both strains used in the 
## mixes (one normal and one hybrid)
bg_files_sub = c("SuSL.Laeth.Laeth.Mix1_LEM2357.mapq30.removedups.proper_paired.deeptools.bin10000.cured.bedgraph",
                 "SuSL.Laeth.Laeth.Mix1_LEM3498.mapq30.removedups.proper_paired.deeptools.bin10000.cured.bedgraph")

## separate code for Laeth samples - myrthe
bg_files_sub = bg_files[4:21]
meta_data_file = '/Users/pmonsieurs/programming/leishmania_susl/data/sample_submission_form_SuSL.20230116.xlsx'
meta_data = read.xlsx(meta_data_file)
meta_data$GS_ID = gsub("105328-001-", "", meta_data$GS_ID)
meta_data

p_combined = plot_layout(ncol = 5)
plot_count = 0
for (bg_file in bg_files_sub) {
  ## extract the sample name
  # sample_name = gsub(".mapq30.removedups.proper_paired.deeptools.bin10000.cured.bedgraph", "", bg_file)
  sample_name = gsub(".mapq30.removedups.proper_paired.deeptools.bin2000.bedgraph", "", bg_file)
  sample_name = gsub("SuSL.Laeth.Laeth.", "", sample_name)
  print(sample_name)
  
  ## only to be performed for specific subset of data
  sample_name = meta_data[meta_data$GS_ID == sample_name,]$Customer.ID
  print(sample_name)
  
  
  plot_count = plot_count + 1
  
  ## set the CPM values
  cpm = read.csv(bg_file, header=FALSE, sep="\t")
  colnames(cpm) = c('chrom', 'start', 'end', 'cpm')
  cpm$sample = sample_name
  
  
  # create empty datafram with one row per chromsosome for
  # which we save the somy value. This somy value is the 
  # median somy value for each of the somy calculations. This
  # basically corresponds to the horizontal line in the boxplots
  # only take the first 36 chromosomes, some contain alos scaffolds
  # which are unfinished. Also take only the CPM values only for those
  # measurements within the first 36 chromosomes
  # change to row 2 to 37 as there is an additional small contig 
  # starting with 000
  somy_df = as.data.frame(unique(cpm$chr)[1:36])
  colnames(somy_df) = c('chrom')
  cpm = cpm[cpm$chrom %in% unique(cpm$chr)[1:36],]
  
  ## plot raw coverage values ##
  # cpm_sub = cpm[,c('chr', 'start', 'end', sample_name)]  
  # colnames(cpm_sub)[4] = 'cpm'
  plot_title = paste0('CPM per chrom - ', sample_name)
  # make_boxplot(cpm_sub, 'cpm', 750, plot_title, FALSE)
  
  ## raw somy values ##
  
  # Raw somy values: add the somy value per window by dividing 
  # each coverage by the median value of the coverage over 
  # all chromosomes
  median_cov = get_median_overall(cpm)
  if (median_cov == 0) {
    print(paste0("median coverage = 0 --> too low?  [", sample_name, "]"))
    next
  }
  
  cpm$raw_somy = 2*cpm[,'cpm']/median_cov
  # plot_title = paste0(sample_name)
  # make_boxplot(cpm_sub, 'raw_somy', 5, plot_title)
  somy_df$somy_raw = calculate_somy(cpm,'raw_somy')
  # asd_raw = calculate_average_somy_deviation(somy_df$somy_raw, somy_ref[, strain])
  # sdc_raw = calculate_somy_difference_count(somy_df$somy_raw, somy_ref[,strain])
  
  
  ## local correction ##
  
  # calculate the chromosome based on the flanking regions. 
  
  flanking = 7
  # somy_values = c()
  # somy_col_name = 'cov'
  chroms =  unique(cpm$chr)[1:36] ## see above: starting from 2...
  cpm$somy_local = NA
  somy_col_name = 'cpm'
  somy_values = c()
  for (i in 1:length(chroms)) {
    chrom = chroms[i]
    flanking_chroms = chroms[max(0,i-flanking):min(length(chroms), i+flanking)]
    cpm_sub= cpm[cpm$chr %in% flanking_chroms,]
    cov_sub_median = get_median_overall(cpm_sub)
    cov_chrom = median(cpm_sub[cpm_sub$chr==chrom, somy_col_name])
    # print(cov_chrom)
    # print(cov_sub_median)
    cov_chrom = cpm_sub[cpm_sub$chr == chrom,'cpm']
    somy_chrom_local =  2*cov_chrom/cov_sub_median
    cpm[cpm$chr == chrom,]$somy_local = somy_chrom_local
    # print(median(somy_chrom_local))
    somy_values = c(somy_values, median(somy_chrom_local))
  }
  # print(somy_values)
  
  ## do some editing on the names to make the plot more readable
  cpm$chrom = gsub("_v4_pilon", "", cpm$chrom)
  cpm$chrom = gsub("LaeL147_", "", cpm$chrom)
  
  # cpm$sample = gsub("SuSL.Lbraz.Lbraz.", "", cpm$sample)
  plot_title = paste0(sample_name)
  
  ## create boxplot and merge the boxplot with the overall patchwork 
  ## object to make summary plot
  p = make_boxplot(cpm, 'somy_local', somy_truth, 5, plot_title)
  if (plot_count == 1) {
    p_combined = p
  }else{
    p_combined = p_combined + p
  }
  
  ## calculate the statistics on the variation on the somy values
  somy_df$somy_local = somy_values
  
}

p_combined


#### Ltrop code ####

## select L tropica from the complete directory. You might want to select 
## only the clinical samples obtained via SureSelect... This can be 
## performed by selecting for files with ".F"
bg_files_sub = bg_files[grep("SuSL.Ltrop", bg_files)]
bg_files_sub = bg_files_sub[grep("\\.F", bg_files_sub)] ## select clinical SuSL samples

## select L tropica from the clinical directory [see settings at parameter]
## these are !!NOT!! sureselect samples but WGS data obtained via
## Othmane and Hasnaa
# bg_files_sub = bg_files

p_combined = plot_layout(ncol = 5)
plot_count = 0
for (bg_file in bg_files_sub) {
  ## extract the sample name
  sample_name = gsub(".mapq30.removedups.proper_paired.deeptools.bin10000.cured.bedgraph", "", bg_file)
  sample_name = gsub("SuSL.Ltrop.Laeth.", "", sample_name)
  print(sample_name)
  plot_count = plot_count + 1
  
  ## set the CPM values
  cpm = read.csv(bg_file, header=FALSE, sep="\t")
  colnames(cpm) = c('chrom', 'start', 'end', 'cpm')
  cpm$sample = sample_name
  print(head(cpm))
  
  
  # create empty datafram with one row per chromsosome for
  # which we save the somy value. This somy value is the 
  # median somy value for each of the somy calculations. This
  # basically corresponds to the horizontal line in the boxplots
  # only take the first 36 chromosomes, some contain alos scaffolds
  # which are unfinished. Also take only the CPM values only for those
  # measurements within the first 36 chromosomes
  # change to row 2 to 37 as there is an additional small contig 
  # starting with 000
  somy_df = as.data.frame(unique(cpm$chr)[1:36])
  colnames(somy_df) = c('chrom')
  cpm = cpm[cpm$chrom %in% unique(cpm$chr)[1:36],]
  
  ## plot raw coverage values ##
  # cpm_sub = cpm[,c('chr', 'start', 'end', sample_name)]  
  # colnames(cpm_sub)[4] = 'cpm'
  plot_title = paste0('CPM per chrom - ', sample_name)
  # make_boxplot(cpm_sub, 'cpm', 750, plot_title, FALSE)
  
  ## raw somy values ##
  
  # Raw somy values: add the somy value per window by dividing 
  # each coverage by the median value of the coverage over 
  # all chromosomes
  c = get_median_overall(cpm)
  if (median_cov == 0) {
    next
  }
  
  cpm$raw_somy = 2*cpm[,'cpm']/median_cov
  # plot_title = paste0(sample_name)
  # make_boxplot(cpm_sub, 'raw_somy', 5, plot_title)
  somy_df$somy_raw = calculate_somy(cpm,'raw_somy')
  # asd_raw = calculate_average_somy_deviation(somy_df$somy_raw, somy_ref[, strain])
  # sdc_raw = calculate_somy_difference_count(somy_df$somy_raw, somy_ref[,strain])
  
  
  ## local correction ##
  
  # calculate the chromosome based on the flanking regions. 
  
  flanking = 7
  # somy_values = c()
  # somy_col_name = 'cov'
  chroms =  unique(cpm$chr)[1:36] ## see above: starting from 2...
  cpm$somy_local = NA
  somy_col_name = 'cpm'
  somy_values = c()
  for (i in 1:length(chroms)) {
    chrom = chroms[i]
    flanking_chroms = chroms[max(0,i-flanking):min(length(chroms), i+flanking)]
    cpm_sub= cpm[cpm$chr %in% flanking_chroms,]
    cov_sub_median = get_median_overall(cpm_sub)
    cov_chrom = median(cpm_sub[cpm_sub$chr==chrom, somy_col_name])
    # print(cov_chrom)
    # print(cov_sub_median)
    cov_chrom = cpm_sub[cpm_sub$chr == chrom,'cpm']
    somy_chrom_local =  2*cov_chrom/cov_sub_median
    cpm[cpm$chr == chrom,]$somy_local = somy_chrom_local
    # print(median(somy_chrom_local))
    somy_values = c(somy_values, median(somy_chrom_local))
  }
  # print(somy_values)
  
  ## do some editing on the names to make the plot more readable
  cpm$chrom = gsub("_v4_pilon", "", cpm$chrom)
  # cpm$sample = gsub("SuSL.Lbraz.Lbraz.", "", cpm$sample)
  plot_title = paste0(sample_name)
  print(plot_title)
  
  ## create boxplot and merge the boxplot with the overall patchwork 
  ## object to make summary plot
  p = make_boxplot(cpm, 'somy_local', somy_truth, 5, plot_title)
  if (plot_count == 1) {
    p_combined = p
  }else{
    p_combined = p_combined + p
  }
  
  ## calculate the statistics on the variation on the somy values
  somy_df$somy_local = somy_values
  
}

p_combined

#### Ltrop code -  manhattan plot ####

bg_files_sub

## normal Ld31 example
bg_file = 'SuSL.Ldon.Ldono_improved.Mix1.mapq30.removedups.proper_paired.deeptools.bin10000.cured.bedgraph'
prev_chrom = "Ld01"

## tropica example
prev_chrom = "LtrL590_01"
bg_file = 'SuSL.Ltrop.Laeth.FJ20_02.mapq30.removedups.proper_paired.deeptools.bin10000.cured.bedgraph'


sample_name = gsub(".mapq30.removedups.proper_paired.deeptools.bin10000.cured.bedgraph", "", bg_file)
sample_name = gsub("SuSL.Ltrop.Laeth.", "", sample_name)
print(sample_name)
plot_count = plot_count + 1

## set the CPM values
cpm = read.csv(bg_file, header=FALSE, sep="\t")
colnames(cpm) = c('chrom', 'start', 'end', 'cpm')
cpm$sample = sample_name
cpm = cpm[cpm$chrom %in% unique(cpm$chr)[1:36],]


bin_width=10000
cpm$cum_pos = seq(1,dim(cpm)[1]*bin_width, bin_width)

cpm$cum_pos = NA
cpm$cum_pos[1] = 1



## intialization
start_chrom = 1

## empty vector
tick_positions = c()
tick_labels = c()


for (i in 2:nrow(cpm)) {
  cpm$cum_pos[i] =  cpm$cum_pos[i-1] + bin_width
  if (cpm$chrom[i] != prev_chrom) {
    print("new chrom!")
    end_chrom = cpm$cum_pos[i]
    tick_pos = (start_chrom + end_chrom)/2
    tick_labels = c(tick_labels, prev_chrom)
    tick_positions = c(tick_positions, tick_pos)
    
    ## new inti
    prev_chrom = cpm$chrom[i]
    start_chrom = cpm$cum_pos[i]
    
  }
}

## add the last label
end_chrom = cpm$cum_pos[i]
tick_pos = (start_chrom + end_chrom)/2
tick_labels = c(tick_labels, prev_chrom)
tick_positions = c(tick_positions, tick_pos)







# plot of the normalized coverage
overall_median = median(cpm$cpm, na.rm=TRUE)
# cpm$somy = 2*(cpm$cov/overall_median)

p2 = ggplot(cpm, aes(x=cum_pos, y=cpm))
p2 = p2 + geom_point(aes(colour=chrom), size=0.50, alpha=0.50)
p2= p2 + scale_color_manual(values = rep(c("#276FBF", "#183059"),
                                         nrow(cpm)))
p2 = p2 + theme(legend.position = "none") 
p2 = p2 + coord_cartesian(ylim=c(0,800))
p2 = p2 + theme_bw()
p2 = p2 + scale_x_continuous(breaks=tick_positions,
                             labels=as.character(tick_labels))
p2 = p2 + theme(axis.text.x = element_text(angle = 90, 
                                           hjust = -.5, 
                                           size=rel(1)))
p2 = p2 + theme(legend.position = "none")

p2




## focus only on chrom 31
cpm=cpm[cpm$chrom == 'LtrL590_31',]
cpm=cpm[cpm$chrom == 'Ld31',]



p2 = ggplot(cpm, aes(x=start, y=cpm))
p2 = p2 + geom_point(aes(colour=chrom), size=1, alpha=0.80)
p2= p2 + scale_color_manual(values = rep(c("#276FBF", "#183059"),
                                         nrow(cpm)))
p2 = p2 + theme(legend.position = "none") 
p2 = p2 + coord_cartesian(ylim=c(0,800))
p2 = p2 + theme_bw()
p2 = p2 + theme(axis.text.x = element_text(angle = 90, 
                                           hjust = -.5, 
                                           size=rel(1)))
p2 = p2 + theme(legend.position = "none")
p2 = p2 + stat_smooth(method="lm", color="red", span=2000, alpha=.50)


p2



##### Additional Code #####


### GC-percentage of chromosome versus somy value. 
## GC percentage of chromosomes is obtained from data from the 
## single cell manuscript (picoplex)
GC_per_chrom_file = '/Users/pmonsieurs/programming/leishmania_singlecell/data/genomes/Leishmania_donovani_16Nov2015beta.gc.csv'
gc_data = read.csv(GC_per_chrom_file, header=TRUE)
gc_vs_somy_data = as.data.frame(gc_data[-37, c('chrom', 'gc_perc')])
gc_vs_somy_data
gc_vs_somy_data = cbind(gc_vs_somy_data, somy_df$somy_local)
colnames(gc_vs_somy_data)[3] = 'somy'
rownames(gc_vs_somy_data) = gc_vs_somy_data$chrom

# remove Ld31
gc_vs_somy_data = gc_vs_somy_data[-which(gc_vs_somy_data$chrom=="Ld31"),]

plot_title = paste0(strain, ' - ', method, ' - somy versus %GC')


p = ggplot(gc_vs_somy_data, aes(x=gc_perc, y=somy))
p = p + geom_point(color="blue", size=2, alpha=0.75)
p = p + stat_smooth(method="lm", color="red", span=2, alpha=.50)
p = p + theme_bw()
p = p + geom_text(aes(label=chrom), hjust=-0.25, vjust=-0.10, size=3)
p = p + coord_cartesian(ylim=c(1.5,2.5))
p = p + ggtitle(plot_title)
p 




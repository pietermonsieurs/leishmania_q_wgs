library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(reshape2)
library(stringr)
library(openxlsx)
library(EnhancedVolcano)
library(VennDiagram)
library(dplyr)
library(tidyr)
library(clusterProfiler)


#### Section 1: Input Data  ####-------------------------

# ---- 1.1: set file names ----
# read in the file names from the directory based on a pattern 
# match. Only the ones -----ending with .csv should be kept. 
src_dir = '/Users/pmonsieurs/programming/leishmania_10X_Q/results/bulk/bwa_htseq/'
out_dir = '/Users/pmonsieurs/programming/leishmania_10X_Q/results/bulk/bwa_htseq/'
gff_file = '/Users/pmonsieurs/programming/leishmania_10X_Q/data/refgenome/Leishmania_donovani_16Nov2015beta.gff3'
gene_id_file = '/Users/pmonsieurs/programming/leishmania_10X/data/blast/Mapping_Sanger_vs_TriTrypDB.csv'
meta_data_file = '/Users/pmonsieurs/programming/leishmania_10X_Q/data/metadata_bulk_RNAseq.xlsx'
setwd(src_dir)

# select the ouput files of htseq count
sample_files <- list.files( path=src_dir, pattern = "*.htseq.strand_no.csv$", full.names = TRUE )
sample_files


# read in the data, and store them in a list, one entry 
# per sample file. The column number is an important parameters
# as it reflects the overall count, sense count or antisense
# count. 
counts.files <- lapply(sample_files, read.table, skip = 0 )
col_number = 2
counts <- as.data.frame( sapply( counts.files, function(x) x[ , col_number ] ) )
head(counts)



## set rownames and column names of the count matrix. For the column names, rely 
## on the meta data file
row.names(counts) <- counts.files[[1]]$V1
meta_data = read.xlsx(meta_data_file)
colnames(meta_data)[6] = 'condition'
meta_data$sample_name = paste0(meta_data$BPK, "_", gsub(" ", "_", meta_data$condition), "_", meta_data$replicate)
meta_data$sample_name
colnames(counts) = meta_data$sample_name
head(counts)

out_file_xlsx = paste0(out_dir, 'raw_read_counts_per_sample.xlsx')
write.xlsx(counts, file = out_file_xlsx, rowNames = TRUE)


## !!!!! remove D4 samples !!!!!!! - only for PCA plots
# counts = counts[,grep("D4", colnames(counts), invert=TRUE)]

# ---- 1.2: create sampleTable ----
# define the sample name and the biological condition and use this 
# information to build a data frame as input for the DESeq2 object
coldata = as.data.frame(colnames(counts))
conditions = paste0(meta_data$BPK, "_", meta_data$condition)
conditions = gsub(" ", "_", conditions)
# conditions = conditions[grep("D4", conditions, invert = TRUE)]
coldata$condition = as.factor(conditions)
colnames(coldata) = c('sample', 'condition')
coldata


# ---- 1.3: create DESeq object with design ----
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design= ~ condition)
dds = DESeq(dds)


# the following selection is done to reduce the number of genes
# with a low expression level. Has *no* effect in the case of quiescence
# experiment as no genes have less than 100 reads over all samples. 
keep <- rowSums(counts(dds)) >= 100
sum(keep)
dds <- dds[keep,]



# ---- 1.4: read GFF file  ----

gff = read.csv(gff_file, sep="\t", comment.char="#", header=FALSE)
unique(gff$V3)
colnames(gff)

## extract the function of a gene
# gff_function = gff[gff$V3 == "protein_coding_gene",]
gff_function = gff
gff_function$description = ""
gff_function$name = ""
gff_function$geneID = ""

for (i in 1:nrow(gff_function)) {
  field9 = str_split_fixed(gff_function[i,9], ";",n=Inf)
  gff_function$gene_id[i] = gsub("ID=", "", field9[1])
  gene_type = str_split_fixed(field9[4], ":", n=Inf)[2]
  if (is.na(gene_type)) {
    gene_type = str_split_fixed(field9[3], ":", n=Inf)[2]
  }
  
  gff_function$gene_type[i] = gene_type
  
  for (j in 1:length(field9)) {
    print(field9[j])
    if(grepl("product=", field9[j])) {
      print(field9)
      description = gsub("product=", "", field9[j])
      description = gsub("term=", "", description)
      
      print(description)
      gff_function$description[i] = description
    }
    
    if(grepl("Name=", field9[j])) {
      name = gsub("Name=", "", field9[j])
      gff_function$name[i] = name
    }
  }
  
}
head(gff_function)
unique(gff_function$gene_type)

# gff = gff[gff$V3 == "exon",]
field9 = str_split_fixed(gff[,9], ";",n=Inf)
gene_id = gsub("ID=", "", field9[,1])

gff = gff[,c(1,3,4,5,7)]
colnames(gff) = c("chrom", "feature", "start", "stop", "strand")
gff$length = as.integer(gff$stop) - as.integer(gff$start) + 1

gff$gene_id = gene_id
head(gff)

index = match(rownames(dds), gff$gene_id)
mcols(dds)$basepairs = gff$length[index]


index = match(rownames(dds), gff$gene_id)
mcols(dds)$basepairs = gff$length[index]
dds_fpkm = fpkm(dds)

# ---- 1.5 read gene IDs -----
gene_data = read.table(gene_id_file)
colnames(gene_data) = c('id_v1', 'id_v2')
gene_data$id_v1 <- sub("\\.1$", "", gene_data$id_v1)



#### Section 2: Quality Control #### 

# ---- 2.1 conversion for QC ---- 
# first make some conversion that make it easier to do some 
# quality control plots
ntd <- normTransform(dds)
vsd <- vst(dds, blind=TRUE)
rld <- rlog(dds, blind=TRUE)
head(assay(vsd), 3)

## save the vsd object to be used as input for the mfuzz clustering
vsd_file = paste0(out_dir, 'vsd_object_for_mfuzz.R')
saveRDS(vsd, file = vsd_file)


# ---- 2.2 PCA / correlation matrix ----
sampleDists <- dist(t(assay(vsd)))

# correlation plot shown as heatmap
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vsd))
colnames(sampleDistMatrix) <- NULL
#colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
p = pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             fontsize=12) #,
#         col=colors)
p
png_file_cor = paste0(out_dir, 'correlation_replicates.png')
ggsave(png_file_cor, p)


## make a PCA plot with all samples. First extract all the PCA plotting
## data as input for ggplot2
pca_df = plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar = round(100 * attr(pca_df, "percentVar"))

## rename the conditions to make them compatible with naming in the 
## manuscript
pca_df$group = gsub("_Pro_LOG_post_PAT", " - Post PAT", pca_df$condition)
pca_df$group = gsub("_Pro_LOG", " - No PAT D2", pca_df$group)
pca_df$group = gsub("_Pro_PAT_D2", " - PAT D2", pca_df$group)
pca_df$group = gsub("_Pro_PAT_D4", " - PAT D4", pca_df$group)

p = ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = group)
) +
  geom_point(size = 3.5, alpha = 0.7) +
  theme_bw() +
  labs(
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)")
  ) +
  scale_x_continuous(expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = 0)) +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    panel.grid = element_blank()
  )

# compute ranges
xr = range(pca_df$PC1)
yr = range(pca_df$PC2)

# add 3% padding
xpad = diff(xr) * 0.03
ypad = diff(yr) * 0.03

p = p +
  coord_cartesian(
    xlim = c(xr[1] - xpad, xr[2] + xpad),
    ylim = c(yr[1] - ypad, yr[2] + ypad)
  )

p

png_file_pca = paste0(out_dir, 'PCA_replicates.png')
ggsave(file = png_file_pca, p,
       dpi=300,
       width=7,
       height=5)




# ---- 2.3 MA plot ----
DESeq2::plotMA(res)





# ---- 2.4: density plots (including FPKM) ----
# density plots of the counts (raw and normalized)
density_raw = melt(log10(counts(dds)))
colnames(density_raw) = c('genename', 'condition', 'count')
p = ggplot(density_raw, aes(x=count, color=condition))
# p = p + geom_density()
p = p + stat_density(geom="line",position="identity", size=1)
p = p + xlim(c(0,5))
p = p + ggtitle("Raw read count per gene")
p = p + labs(x="log 10 read count per gene")
p = p + theme_bw()
p = p + theme(axis.text=element_text(size=12),
              axis.title=element_text(size=16,face="bold"),
              title=element_text(size=16))
p

dens_plot_raw_file = paste0(out_dir, 'density_plot_raw.png')
ggsave(dens_plot_raw_file, p)

density_raw = melt(log10(counts(dds, normalized=TRUE)))
colnames(density_raw) = c('genename', 'condition', 'count')
p = ggplot(density_raw, aes(x=count, color=condition))
p = p + stat_density(geom="line",position="identity", size=1)
p = p + xlim(c(0,5))
p = p + theme_bw()
p = p + ggtitle("Normalized read count per gene")
p = p + labs(x="log 10 read count per gene")
p = p + theme(axis.text=element_text(size=12),
              axis.title=element_text(size=16,face="bold"),
              title=element_text(size=16))
p

dens_plot_norm_file = paste0(out_dir, 'density_plot_normalized.png')
ggsave(dens_plot_norm_file, p)

# FPKM values and density values
dds_fpkm = fpkm(dds)
fpkm_melted = melt(dds_fpkm)
colnames(fpkm_melted) = c('genename', 'condition', 'fpkm')

p = ggplot(fpkm_melted, aes(x=fpkm, color=condition))
p = p + geom_density()
p = p + xlim(c(0,5))
p = p + theme_bw()
p

dens_plot_norm_file = paste0(out_dir, 'density_plot_fpkm.png')
ggsave(dens_plot_norm_file, p)



#### Section 3: Differential expression  ####

# ---- 3.1 compare Pro LOG vs Pro PAT D2 ---- 

res = results(dds, contrast = c("condition", "BPK275_Pro_PAT_D2", "BPK275_Pro_LOG"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = paste0("BPK275_Pro_PAT_D2_VS_Pro_LOG", c('_fc', '_pval', '_fdr'))
df_out = res
head(df_out)

res = results(dds, contrast = c("condition", "BPK026_Pro_PAT_D2", "BPK026_Pro_LOG"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = paste0("BPK026_Pro_PAT_D2_VS_Pro_LOG", c('_fc', '_pval', '_fdr'))
df_out = cbind.data.frame(df_out,res)
head(df_out)

# ---- 3.2 compare Pro PAT D4 vs Pro PAT D2 for BPK275 ---- 
res = results(dds, contrast = c("condition", "BPK275_Pro_PAT_D4", "BPK275_Pro_PAT_D2"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = paste0("BPK275_Pro_PAT_D4_VS_Pro_PAT_D2", c('_fc', '_pval', '_fdr'))
df_out = cbind.data.frame(df_out,res)
head(df_out)

## additional request: compare Pro PAT D4 vs Pro PAT D2 for BPK275 but only for 
## the two high quality D2 samples, those are stored in column 7 instead of 6 
## in the metadata file
meta_data2 = read.xlsx(meta_data_file)
colnames(meta_data2)[7] = 'condition'
meta_data2$sample_name = paste0(meta_data2$BPK, "_", gsub(" ", "_", meta_data2$condition), "_", meta_data2$replicate)
meta_data2$sample_name
counts2 = counts
colnames(counts2) = meta_data2$sample_name
head(counts2)

# define the sample name and the biological condition and use this 
# information to build a data frame as input for the DESeq2 object
coldata2 = as.data.frame(colnames(counts2))
conditions2 = paste0(meta_data2$BPK, "_", meta_data2$condition)
conditions2 = gsub(" ", "_", conditions2)
coldata2$condition = as.factor(conditions2)
colnames(coldata2) = c('sample', 'condition')
coldata2


dds2 <- DESeqDataSetFromMatrix(countData = counts2,
                              colData = coldata2,
                              design= ~ condition)
keep <- rowSums(counts(dds2)) >= 100
dds2 <- dds2[keep,]
dds2 = DESeq(dds2)

res = results(dds2, contrast = c("condition", "BPK275_Pro_PAT_D4", "BPK275_Pro_PAT_D2"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = paste0("BPK275_Pro_PAT_D4_VS_Pro_PAT_D2_reduced", c('_fc', '_pval', '_fdr'))
df_out = cbind.data.frame(df_out,res)
head(df_out)

## additinonal request; compare BPK275 PRO LOG and BPK275 PRO PAT D4
res = results(dds, contrast = c("condition", "BPK275_Pro_PAT_D4", "BPK275_Pro_LOG"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = paste0("BPK275_Pro_PAT_D4_VS_Pro_LOG", c('_fc', '_pval', '_fdr'))
df_out = cbind.data.frame(df_out,res)
head(df_out)



# ---- 3.3 compare Pro LOG post PAT vs Pro LOG ---- 

res = results(dds, contrast = c("condition", "BPK275_Pro_LOG_post_PAT", "BPK275_Pro_LOG"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = paste0("BPK275_Pro_LOG_post_PAT_VS_Pro_LOG", c('_fc', '_pval', '_fdr'))
df_out = cbind.data.frame(df_out,res)
head(df_out)

res = results(dds, contrast = c("condition", "BPK026_Pro_LOG_post_PAT", "BPK026_Pro_LOG"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = paste0("BPK026_Pro_LOG_post_PAT_VS_Pro_LOG", c('_fc', '_pval', '_fdr'))
df_out = cbind.data.frame(df_out,res)
head(df_out)


# ---- 3.4 compare between both strains ---- 
res = results(dds, contrast = c("condition", "BPK275_Pro_LOG", "BPK026_Pro_LOG"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = paste0("Pro_LOG_BPK275_vs_BPK026", c('_fc', '_pval', '_fdr'))
df_out = cbind.data.frame(df_out,res)
head(df_out)


res = results(dds, contrast = c("condition", "BPK275_Pro_PAT_D2", "BPK026_Pro_PAT_D2"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = paste0("Pro_PAT_D2_BPK275_vs_BPK026", c('_fc', '_pval', '_fdr'))
df_out = cbind.data.frame(df_out,res)
head(df_out)


res = results(dds, contrast = c("condition", "BPK275_Pro_LOG_post_PAT", "BPK026_Pro_LOG_post_PAT"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = paste0("Pro_LOG_post_PAT_BPK275_vs_BPK026", c('_fc', '_pval', '_fdr'))
df_out = cbind.data.frame(df_out,res)
head(df_out)



# ---- 3.5 add additional columsn (descript + FPKM) ----

df_out$tritrypDB_id = gene_data[match(rownames(df_out),gene_data$id_v1),]$id_v2
df_out$name = gff_function[match(rownames(df_out), gff_function$gene_id),]$name
df_out$description = gff_function[match(rownames(df_out), gff_function$gene_id),]$description
df_out$gene_type = gff_function[match(rownames(df_out), gff_function$gene_id),]$gene_type

df_out$chrom = gff_function[match(rownames(df_out), gff_function$gene_id),1]
df_out$start = gff_function[match(rownames(df_out), gff_function$gene_id),4]
df_out$stop =  gff_function[match(rownames(df_out), gff_function$gene_id),5]
df_out$strand =  gff_function[match(rownames(df_out), gff_function$gene_id),7]

## add the FPKM value
dim(dds_fpkm)
dim(df_out)
sum(rownames(dds_fpkm) == rownames(df_out))

df_out = cbind.data.frame(df_out, dds_fpkm)
head(df_out)


excel_file = paste0(out_dir, 'DESeq2.xlsx')
write.xlsx(df_out, file=excel_file, rowNames=TRUE)

## save output also to a RDS file
rds_file = paste0(out_dir, 'DESeq2_df_out.rds')
saveRDS(df_out, file = rds_file)
# df_out = readRDS(rds_file)

## plot the distribution of the RPKM values to see what are reasonable ranges
## and what is higher or low expressed
rpkm_df = df_out[,grep("_R", colnames(df_out))]
rpkm_plot = melt(rpkm_df)
rpkm_plot$strain <- sub("_.*", "", rpkm_plot$variable)
rpkm_plot$condition <- sub("^(BPK\\d+_)?(.*?)(_R\\d+)?$", "\\2", rpkm_plot$variable)

head(rpkm_plot)

ggplot(data = rpkm_plot, aes(x=value)) + 
  geom_density(aes(colour = strain), n=5*4096) + 
  facet_wrap( ~ condition) + 
  coord_cartesian(xlim = c(0,100)) + 
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  theme_bw()

ggplot(data = rpkm_plot, aes(x=value)) + 
  geom_density(aes(colour = strain, linetype = condition), n=5*4096) + 
  coord_cartesian(xlim = c(0,100)) + 
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  theme_bw()


## extract some statistics on up and downregulation
results_list <- list()

for (comparison in grep("_fc", colnames(df_out), value = TRUE)) {
  pval_column <- gsub("_fc", "_fdr", comparison)
  row_res <- c(comparison = comparison)
  
  for (cutoff in c(0.58, 0.75, 1)) {
    up <- nrow(df_out[df_out[, comparison] > cutoff & df_out[, pval_column] < 0.05, , drop = FALSE])
    down <- nrow(df_out[df_out[, comparison] < -cutoff & df_out[, pval_column] < 0.05, , drop = FALSE])
    
    row_res[paste0("up_", cutoff)] <- up
    row_res[paste0("down_", cutoff)] <- down
  }
  
  results_list[[comparison]] <- row_res
}

## combine into a single data frame
results_df <- do.call(rbind, results_list)
results_df <- as.data.frame(results_df, stringsAsFactors = FALSE)

out_file = paste0(out_dir, 'stats_up_and_down_regulated.xlsx')
write.xlsx(results_df, file = out_file)



# ---- 3.6 Volcano plots -----

volcano_strain = 'BPK026' 
volcano_strain = 'BPK275' 

EnhancedVolcano(df_out,
                lab = rownames(res),                   # gene labels
                x = paste0(volcano_strain, '_Pro_PAT_D2_VS_Pro_LOG_fc'),
                y = paste0(volcano_strain, '_Pro_LOG_post_PAT_VS_Pro_LOG_fdr'),
                pCutoff = 0.05,                        # p-value cutoff
                FCcutoff = 1.0,                        # fold change cutoff
                pointSize = 1.0,
                labSize = 3.0,
                title = paste0("Volcano plot:", volcano_strain, " PAT D2 versus Pro LOG(ref)"),
                subtitle = "DESeq2 results",
                caption = "Thresholds: abs(log2FC) > 1 and p < 0.05",
                legendLabels = c("NS","Log2FC","p-value","p-value & Log2FC"),
                legendPosition = 'right',
                colAlpha = 0.9,
                drawConnectors = TRUE,
                widthConnectors = 0.5) + coord_cartesian(ylim = c(0,10))

# ---- 3.7 GO enrichment analysis -----

## create a GO library based on the GO classes as obtained from the TriTrypDB
## for the TriTrypDB genome. The gene names still needed to be converted from 
## TriTrypDB to the v2 Sanger genome, using the script bulk_GO_convert_gaf.py
## That is also why the "v2" is added to the gaf file 
gaf_file = '/Users/pmonsieurs/programming/leishmania_10X_Q/data/refgenome/TriTrypDB-68_LdonovaniBPK282A1v2_GO.gaf'
gaf_file = '/Users/pmonsieurs/programming/leishmania_10X_Q/data/refgenome/TriTrypDB-68_LdonovaniBPK282A1v2_Curated_GO.gaf'
go_obo_file = '/Users/pmonsieurs/programming/plasmodium_hiveseq/data/refgenome/go-basic.obo'

gaf = read.delim(gaf_file, 
                 header = FALSE, 
                 comment.char = "!", 
                 stringsAsFactors = FALSE)
head(gaf)

## create a TERM2GENE dataframe that can be used in the enricher function
term2gene <- gaf[, c(5, 2)]
colnames(term2gene) <- c("GO", "gene_id")
term2gene <- unique(term2gene)
head(term2gene)

## create an additional data frame to have more explanation on the GO term - 
## otherwise no useful information is printed. The GO information is coming from
## a generic obo file. Extract the ids and the names, and they should be in the 
## same order in the obo file. First do a grep, then replace
obo <- readLines(go_obo_file)
ids <- grep("^id: GO:", obo, value = TRUE)
names <- grep("^name: ", obo, value = TRUE)
go_ids <- sub("^id: ", "", ids)
go_names <- sub("^name: ", "", names[seq_along(go_ids)]) 
term2name <- data.frame(GO = go_ids, Description = go_names, stringsAsFactors = FALSE)


## !! test example !! ##
gene_list = rownames(df_out)[c(100,1000,4000,4005,5000)]
ego <- enricher(
  gene = gene_list,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pAdjustMethod = "none", # "BH"
  qvalueCutoff = 0.05
)


head(ego)
barplot(ego, showCategory=20)
dotplot(ego, showCategory=20)


## get the list of differentially expressed genes for different combinations
logfc_cutoff = .58
fdr_cutoff = .05

diff_expressed_list = list()

diff_expressed_list[['BPK275_Pro_PAT_D2_all']] = rownames(df_out[abs(df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fc) > logfc_cutoff &
                                                                    df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fdr <= fdr_cutoff,])
diff_expressed_list[['BPK026_Pro_PAT_D2_all']] = rownames(df_out[abs(df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fc) > logfc_cutoff &
                                                                    df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fdr <= fdr_cutoff,])
diff_expressed_list[['intersect_all']] = intersect(diff_expressed_list[['BPK275_Pro_PAT_D2_all']],
                                                    diff_expressed_list[['BPK026_Pro_PAT_D2_all']])

diff_expressed_list[['BPK275_Pro_PAT_D2_up']] = rownames(df_out[df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fc > logfc_cutoff &
                                                      df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fdr <= fdr_cutoff,])
diff_expressed_list[['BPK026_Pro_PAT_D2_up']] = rownames(df_out[df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fc > logfc_cutoff &
                                                               df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fdr <= fdr_cutoff,])
diff_expressed_list[['intersect_up']] = intersect(diff_expressed_list[['BPK275_Pro_PAT_D2_up']],
                                                          diff_expressed_list[['BPK026_Pro_PAT_D2_up']])

diff_expressed_list[['BPK275_Pro_PAT_D2_down']] = rownames(df_out[df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fc < - logfc_cutoff &
                                                                  df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fdr <= fdr_cutoff,])
diff_expressed_list[['BPK026_Pro_PAT_D2_down']] = rownames(df_out[df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fc < - logfc_cutoff &
                                                                  df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fdr <= fdr_cutoff,])
diff_expressed_list[['intersect_down']] = intersect(diff_expressed_list[['BPK275_Pro_PAT_D2_down']],
                                                  diff_expressed_list[['BPK026_Pro_PAT_D2_down']])



## run the enrichment analysis per cluster. 
enrich_list <- lapply(names(diff_expressed_list), function(condition) {
  print(condition)
  gene_vec <- diff_expressed_list[[condition]]
  print(length(gene_vec))
  
  ego <- enricher(
    gene = gene_vec,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pAdjustMethod = "BH", # BH
    qvalueCutoff = 1,
    pvalueCutoff = 0.05,
    minGSSize = 5
  )
  
  if (nrow(ego) > 0) {
    df <- as.data.frame(ego)
    df$cluster <- condition
    return(df)
  } else {
    return(NULL)  # skip clusters with no enriched terms
  }
})

lapply(enrich_list, function(x) x[, 1:4])

## the enrich_list is now still a list, with 9 different items, each item a 
## data frame. Now combine all in one data frame
enrich_df <- do.call(rbind, lapply(enrich_list, as.data.frame))
enrich_df[grepl("utophagy", enrich_df$Description),]


top_terms = enrich_df


## factor for plotting order. This ensures that the ordering is first is based
## on the cluster, and then ordered based on the p-value 
top_terms <- top_terms %>%
  arrange(cluster, p.adjust) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

go_dotplot = ggplot(top_terms, aes(x = cluster, y = Description, size = Count, color = -log10(p.adjust))) +
  geom_point() +
  scale_color_gradient(low = "skyblue", high = "red") +
  theme_bw() +
  labs(
    x = "Cluster",
    y = "GO Term",
    color = "-log10(adj.p)",
    size = "Gene count",
    title = "GO enrichment"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

go_dotplot

go_dotplot_file = paste0(out_dir, 'go_dotplot.png')
ggsave(file = go_dotplot_file,
       plot = go_dotplot,
       width = 9,
       height = 12,
       dpi = 300)

go_enricher_file = paste0(out_dir, 'enricher_go.xlsx')
write.xlsx(enrich_df, file = go_enricher_file)


## focus on the autophagy pathway
go_autophagy = term2name[term2name$Description == "autophagy",]$GO
go_autophagy_genes = term2gene[term2gene$GO == go_autophagy,]$gene_id

df_out_autophagy = df_out[rownames(df_out) %in% go_autophagy_genes,]
plot_data = df_out_autophagy[, grep("_fc", colnames(df_out_autophagy))]
rownames(plot_data) = paste0(rownames(plot_data), "::", gff_function[match(rownames(plot_data), gff$gene_id),]$description)
#plot_mat <- as.matrix(apply(plot_data, 2, as.numeric))

breaks = seq(-2,2,length.out = 100)
pheatmap(plot_data,
         breaks = breaks,
         color = colorRampPalette(c("blue", "white", "red"))(99),
         display_numbers = TRUE,
         number_format = "%.2f"
         )





go_autophagy = term2name[term2name$Description == "autophagy",]$GO
term2gene[term2gene$GO == go_autophagy,]$gene_id



# ---- 3.8 KEGG enrichment analysis -----

## create a KEGG library similarly as was done for the GO analysis, and try
## to reproduce the GO structure of term2gene and term2name
gaf_file = '/Users/pmonsieurs/programming/leishmania_10X_Q/data/refgenome/TriTrypDB-48_LdonovaniBPK282A1.KEGG.gmt'
go_obo_file = '/Users/pmonsieurs/programming/leishmania_10X_Q/data/refgenome/TriTrypDB-48_LdonovaniBPK282A1.KEGG.gmt'
kegg_to_skip_file = '/Users/pmonsieurs/programming/leishmania_10X_Q/results/bulk/bwa_htseq/go/kegg_pathways_to_be_skipped.list'
gaf = read.delim(gaf_file, 
                 header = FALSE, 
                 comment.char = "!", 
                 stringsAsFactors = FALSE)
head(gaf)

## read in kegg names to skip becasue not relevant for Leishmania
kegg_to_skip = read.csv(kegg_to_skip_file)[,1]
kegg_to_skip = paste0("ko", sprintf("%05d", kegg_to_skip))

## defined the maximum number of columsn in the KEGG datafile
counts <- count.fields(gaf_file, sep = "\t")
table(counts)
maxn <- max(counts) ## maximum number of columns present
cc <- rep("character", maxn)

## count fields per line (tab-separated), and select the maximum number to 
## define the number of columns
obo_full <- read.table(
  gaf_file,
  header = FALSE,
  sep = "\t",
  comment.char = "!",
  stringsAsFactors = FALSE,
  fill = TRUE,
  quote = "",
  colClasses = cc
)
dim(obo_full)

## create term2name
term2name <- obo_full[, 1:2]
colnames(term2name) = c("GO", "Description")

## create term2gene
obo_reduced = obo_full[,-2]
obo_melt = melt(obo_reduced, id.vars = "V1")
unique(obo_melt[,2])

obo_melt = obo_melt[obo_melt[,3] != "",]
unique(obo_melt[,3])
dim(obo_melt)

## map TriTrypDB Ids to Sanger IDs
mapping_file = '/Users/pmonsieurs/programming/leishmania_10X_Q/data/refgenome/Mapping_TriTrypDB_vs_Sanger.csv'
mapping_data = read.table(mapping_file)
colnames(mapping_data) = c('TriTrypDB', 'Sanger')
mapping_data$TriTrypDB = gsub("\\.1\\.1", "", mapping_data$TriTrypDB)
mapping_data$Sanger = gsub("\\.1", "", mapping_data$Sanger)
head(mapping_data)

obo_melt$sanger = mapping_data[match(obo_melt[,3], mapping_data$TriTrypDB),]$Sanger
obo_melt = obo_melt[!is.na(obo_melt[,4]),]
dim(obo_melt)
head(obo_melt)

term2gene = obo_melt[,c(1,4)]
colnames(term2gene) <- c("GO", "gene_id")


## remove the kegg pathways which are not relevant
term2gene = term2gene[! term2gene$GO %in% kegg_to_skip,]
term2name = term2name[! term2name$GO %in% kegg_to_skip,]

## !! test example !! ##
gene_list = rownames(df_out)[c(100,1000,4000,4005,5000)]
ego <- enricher(
  gene = gene_list,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pAdjustMethod = "none", # "BH"
  qvalueCutoff = 0.05
)


head(ego)
barplot(ego, showCategory=20)
dotplot(ego, showCategory=20)


## get the list of differentially expressed genes for different combinations
# logfc_cutoff = .75
logfc_cutoff = 1
fdr_cutoff = 0.05

diff_expressed_list = list()

diff_expressed_list[['BPK275_Pro_PAT_D2_all']] = rownames(df_out[abs(df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fc) > logfc_cutoff &
                                                                   df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fdr <= fdr_cutoff,])
diff_expressed_list[['BPK026_Pro_PAT_D2_all']] = rownames(df_out[abs(df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fc) > logfc_cutoff &
                                                                   df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fdr <= fdr_cutoff,])
diff_expressed_list[['intersect_all']] = intersect(diff_expressed_list[['BPK275_Pro_PAT_D2_all']],
                                                   diff_expressed_list[['BPK026_Pro_PAT_D2_all']])

diff_expressed_list[['BPK275_Pro_PAT_D2_up']] = rownames(df_out[df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fc > logfc_cutoff &
                                                                  df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fdr <= fdr_cutoff,])
diff_expressed_list[['BPK026_Pro_PAT_D2_up']] = rownames(df_out[df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fc > logfc_cutoff &
                                                                  df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fdr <= fdr_cutoff,])
diff_expressed_list[['intersect_up']] = intersect(diff_expressed_list[['BPK275_Pro_PAT_D2_up']],
                                                  diff_expressed_list[['BPK026_Pro_PAT_D2_up']])

diff_expressed_list[['BPK275_Pro_PAT_D2_down']] = rownames(df_out[df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fc < - logfc_cutoff &
                                                                    df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fdr <= fdr_cutoff,])
diff_expressed_list[['BPK026_Pro_PAT_D2_down']] = rownames(df_out[df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fc < - logfc_cutoff &
                                                                    df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fdr <= fdr_cutoff,])
diff_expressed_list[['intersect_down']] = intersect(diff_expressed_list[['BPK275_Pro_PAT_D2_down']],
                                                    diff_expressed_list[['BPK026_Pro_PAT_D2_down']])



## run the enrichment analysis per cluster. 
enrich_list <- lapply(names(diff_expressed_list), function(condition) {
  print(condition)
  gene_vec <- diff_expressed_list[[condition]]
  print(length(gene_vec))
  
  ego <- enricher(
    gene = gene_vec,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pAdjustMethod = "BH", # BH
    qvalueCutoff = .5,
    pvalueCutoff = 0.05,
    minGSSize = 2
  )
  
  if (nrow(ego) > 0) {
    df <- as.data.frame(ego)
    df$cluster <- condition
    return(df)
  } else {
    return(NULL)  # skip clusters with no enriched terms
  }
})

lapply(enrich_list, function(x) x[, 1:4])

## the enrich_list is now still a list, with 9 different items, each item a 
## data frame. Now combine all in one data frame
enrich_df <- do.call(rbind, lapply(enrich_list, as.data.frame))
enrich_df[grepl("tophagy", enrich_df$Description),]


top_terms = enrich_df


## factor for plotting order. This ensures that the ordering is first is based
## on the cluster, and then ordered based on the p-value 
top_terms <- top_terms %>%
  arrange(cluster, p.adjust) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

kegg_dotplot = ggplot(top_terms, aes(x = cluster, y = Description, size = Count, color = -log10(p.adjust))) +
  geom_point() +
  scale_color_gradient(low = "skyblue", high = "red") +
  theme_bw() +
  labs(
    x = "Cluster",
    y = "KEGG pathway",
    color = "-log10(adj.p)",
    size = "Gene count",
    title = "KEGG enrichment across clusters"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

kegg_dotplot

kegg_dotplot_file = paste0(out_dir, 'kegg_dotplot.png')
ggsave(file = kegg_dotplot_file,
       plot = kegg_dotplot,
       width = 9,
       height = 12,
       dpi = 300)


kegg_enricher_file = paste0(out_dir, 'enricher_kegg_logfc', logfc_cutoff, '.xlsx')
write.xlsx(enrich_df, file = kegg_enricher_file)



#### 4. custom analysis ####

# ---- 4.1 common up & down at D2 in both strains ----
cutoff_fc = 1 ## also try 0.58 and 1 / 0.75
cutoff_pval = 0.05
up_BPK275 = rownames(df_out[df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fc >= cutoff_fc & df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fdr < cutoff_pval,])
up_BPK026 = rownames(df_out[df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fc >= cutoff_fc & ! is.na(df_out$BPK026_Pro_LOG_post_PAT_VS_Pro_LOG_fdr) & df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fdr < cutoff_pval,])
up_intersect = intersect(up_BPK275, up_BPK026)
length(up_intersect)

down_BPK275 = rownames(df_out[df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fc <= - cutoff_fc & df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fdr < cutoff_pval,])
down_BPK026 = rownames(df_out[df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fc <= - cutoff_fc & ! is.na(df_out$BPK026_Pro_LOG_post_PAT_VS_Pro_LOG_fdr) & df_out$BPK026_Pro_PAT_D2_VS_Pro_LOG_fdr < cutoff_pval,])
down_intersect = intersect(down_BPK275, down_BPK026)
length(down_intersect)

## make venn diagram of overlap between 0.75 up- and down-regulated genes
myCol <- brewer.pal(3, "Pastel2")

## do scaled and non-scaled
scaled = TRUE
if (scaled) {
  scaled_text = "_scaled"
}else{
  scaled_text = "_nonscaled"
}

venn.diagram(
  

  x = list(up_BPK026, up_BPK275),
  filename = paste0(out_dir, 'venn_diagramm_up', scaled_text,'.png'),

  # x = list(down_BPK026, down_BPK275),
  # filename = paste0(out_dir, 'venn_diagramm_down', scaled_text,'.png'),

  
  category.names = c("BPK026" , "BPK275"),
  
  
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  
  # # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],
  scaled = scaled,

  # # Numbers
  cex = .45,
  fontface = "bold",
  fontfamily = "sans",
  
  
  # Set names
  cat.cex = 0.8,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  # cat.pos = c(-25, 18),
  cat.pos = c(0, 0),
  
  cat.dist = c(0.03, 0.03),
  cat.fontfamily = "sans",
  #rotation = 1
)

## save the output for both up and downregulated genes
## save UP 
excel_file = paste0(out_dir, 'DESeq2_intersect_up_cutoff_', cutoff_fc, '.xlsx')
df_up = df_out[rownames(df_out) %in% up_intersect, ]
dim(df_up)
write.xlsx(df_up, file=excel_file, rowNames=TRUE)

## save DOWN
excel_file = paste0(out_dir, 'DESeq2_intersect_down_cutoff_', cutoff_fc, '.xlsx')
df_down = df_out[rownames(df_out) %in% down_intersect, ]
dim(df_down)
write.xlsx(df_down, file=excel_file, rowNames=TRUE)


## save the output for BPK026 separately
strain = 'BPK026'
excel_file = paste0(out_dir, 'DESeq2_',strain,'_up_cutoff_', cutoff_fc, '.xlsx')
df_up = df_out[rownames(df_out) %in% up_BPK026, ]
dim(df_up)
write.xlsx(df_up, file=excel_file, rowNames=TRUE)

## save DOWN
excel_file = paste0(out_dir, 'DESeq2_',strain,'_down_cutoff_', cutoff_fc, '.xlsx')
df_down = df_out[rownames(df_out) %in% down_BPK026, ]
dim(df_down)
write.xlsx(df_down, file=excel_file, rowNames=TRUE)

## save the output for BPK275 separately
strain = 'BPK275'
excel_file = paste0(out_dir, 'DESeq2_',strain,'_up_cutoff_', cutoff_fc, '.xlsx')
df_up = df_out[rownames(df_out) %in% up_BPK275, ]
dim(df_up)
write.xlsx(df_up, file=excel_file, rowNames=TRUE)

## save DOWN
excel_file = paste0(out_dir, 'DESeq2_',strain,'_down_cutoff_', cutoff_fc, '.xlsx')
df_down = df_out[rownames(df_out) %in% down_BPK275, ]
dim(df_down)
write.xlsx(df_down, file=excel_file, rowNames=TRUE)




# ----- 4.2 chromosome 5 fold change ----
df_ld05 = df_out[df_out$chrom == "Ld05",]
df_ld05_plot = df_ld05[, c('chrom', 'start'),]
df_ld05_plot = cbind.data.frame(df_ld05_plot, df_ld05[,grep("fc", colnames(df_out))])
head(df_ld05_plot)

df_ld05_melt = melt(df_ld05_plot)
colnames(df_ld05_melt) = c('chrom', 'start', 'condition_full', 'log2fc')
df_ld05_melt$strain = sub("_.*", "", df_ld05_melt$condition)
df_ld05_melt$condition <- sub("^[^_]*_", "", df_ld05_melt$condition_full)


ggplot(df_ld05_melt, aes(x=start, y=log2fc)) + 
  geom_point(aes(colour = condition)) + 
  theme_bw() + 
  facet_wrap(~ strain, ncol=1) + 
  coord_cartesian(ylim = c(-2,2))
  

mean(df_ld05$BPK275_Pro_LOG_post_PAT_VS_Pro_LOG_fc, na.rm=TRUE)


# ---- 4.3 drug resistance genes ----
dr_file = '/Users/pmonsieurs/programming/leishmania_susl/data/drug_resistance_Ldon_v2.xlsx'
dr_data = read.xlsx(dr_file)
df_dr = df_out[rownames(df_out) %in% dr_data$Sanger_id, ]
df_dr = df_dr[, grep("_R", colnames(df_dr))]

dr_data$print_name = paste0(dr_data$Drug_code, "_", dr_data$Gene_code, "_", dr_data$Sanger_id)
rownames(df_dr) = dr_data[match(rownames(df_dr), dr_data$Sanger_id),]$print_name

pheatmap(df_dr)
pheatmap(df_dr, scale="row")

plot_data = melt(as.matrix(df_dr))
colnames(plot_data) = c('gene', 'sample', 'RPKM')
plot_data$condition = sub("_R[0-9]+$", "", plot_data$sample)
plot_data$strain = sub("_.*", "", plot_data$sample)

## make a factor / level from the condition so you can sort the X-axis to your
## own preference
plot_data$condition <- factor(plot_data$condition,
                              levels = c("BPK026_Pro_LOG", "BPK026_Pro_PAT_D2", "BPK026_Pro_LOG_post_PAT", "BPK275_Pro_LOG", 
                                         "BPK275_Pro_PAT_D2", "BPK275_Pro_PAT_D4", "BPK275_Pro_LOG_post_PAT"))

head(plot_data)

ggplot(plot_data, aes(x=condition, y=RPKM)) + 
  geom_point(aes(color=strain)) + 
  facet_wrap( ~ gene, ncol=4, scale = "free_y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# ---- 4.4 not DE at D2 but DE at D4 ----
cutoff_fc = .75
cutoff_pval = 0.05
not_DE_D2 = rownames(df_out[abs(df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fc) <= cutoff_fc & df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fdr > cutoff_pval,])
up_D4 = rownames(df_out[df_out$BPK275_Pro_PAT_D4_VS_Pro_LOG_fc >= cutoff_fc & df_out$BPK275_Pro_PAT_D4_VS_Pro_LOG_fdr <= cutoff_pval,])
genes_intersect = intersect(not_DE_D2, up_D4)

subset_data = df_out[rownames(df_out) %in% genes_intersect,]
subset_data = subset_data[, grep("_R", colnames(subset_data))]
plot_data = melt(as.matrix(subset_data))
colnames(plot_data) = c('gene', 'replicate', 'rpkm')
plot_data$condition <- sub("_R[0-9]+$", "", plot_data$replicate)
plot_data$strain <- sub("^(BPK[0-9]+).*", "\\1", plot_data$condition)
plot_data$condition_nostrain <- sub("^BPK[0-9]+_", "", plot_data$condition)
plot_data$condition_nostrain = factor(plot_data$condition_nostrain, 
                                      levels= c("Pro_LOG","Pro_PAT_D2",  "Pro_PAT_D4", "Pro_LOG_post_PAT"))


ggplot(plot_data, aes(x=condition_nostrain, y=rpkm)) + 
  geom_point(aes(color=strain), size=2) + 
  facet_wrap(~ gene, scale = 'free_y') + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))




# ---- 4.5 RPKM of ABC transporters ----
abc_transporter_file = '/Users/pmonsieurs/programming/leishmania_10X_Q/data/abctransporters.xlsx'
abc_data = read.xlsx(abc_transporter_file, sheet = 'abc2')
abc_data$Ldon_ortholog_clean = gsub(".1$", "", abc_data$Ldon_ortholog)
head(abc_data)

# count_data = df_out[,24:49]
## only select the columns with RPKM values, they have replicates, so should 
## contain the "_R" pattern
count_data = df_out[,grep("_R", colnames(df_out))]
count_data_abc = count_data[rownames(count_data) %in% abc_data$Ldon_ortholog_clean,]

breaks = seq(1,50,.5)
pheatmap(count_data_abc, 
         breaks=breaks)


## check wheter some of the ABC transporters are up in D4 but no in D2
cutoff_fc = .75
cutoff_pval = 0.05
not_DE_D2 = rownames(df_out[abs(df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fc) <= cutoff_fc & df_out$BPK275_Pro_PAT_D2_VS_Pro_LOG_fdr > cutoff_pval,])
up_D4 = rownames(df_out[df_out$BPK275_Pro_PAT_D4_VS_Pro_LOG_fc >= cutoff_fc & df_out$BPK275_Pro_PAT_D4_VS_Pro_LOG_fdr <= cutoff_pval,])
genes_intersect = intersect(not_DE_D2, up_D4)

intersect(genes_intersect, abc_data$Ldon_ortholog_clean)

## remove the D4 replicates !!!
count_data_abc = count_data_abc[, grep("D4", colnames(count_data_abc), invert = TRUE)]


## plot the RPKM values per ABC transporter for all ABC transporters
# plot_data = melt(as.matrix(count_data_abc[, 7:ncol(count_data_abc)]))
plot_data = melt(as.matrix(count_data_abc))
colnames(plot_data) = c('gene', 'sample', 'RPKM')
plot_data$condition = sub("_R[0-9]+$", "", plot_data$sample)
plot_data$strain = sub("_.*", "", plot_data$sample)

## make a factor / level from the condition so you can sort the X-axis to your
## own preference
plot_data$condition <- factor(plot_data$condition,
                              levels = c("BPK026_Pro_LOG", "BPK026_Pro_PAT_D2", "BPK026_Pro_LOG_post_PAT", "BPK275_Pro_LOG", 
                                         "BPK275_Pro_PAT_D2", "BPK275_Pro_PAT_D4", "BPK275_Pro_LOG_post_PAT"))
plot_data$condition_nostrain <- sub("^BPK[0-9]+_", "", plot_data$condition)
plot_data$condition_nostrain = factor(plot_data$condition_nostrain, 
                                      levels= c("Pro_LOG","Pro_PAT_D2",  "Pro_PAT_D4", "Pro_LOG_post_PAT"))



## do renaming and repeat the full analysis
plot_data$condition_nostrain_clean = gsub("Pro_LOG_post_PAT", "Post PAT", plot_data$condition_nostrain)
plot_data$condition_nostrain_clean = gsub("Pro_LOG", "No PAT D2", plot_data$condition_nostrain_clean)
plot_data$condition_nostrain_clean = gsub("Pro_PAT_D2", "PAT D2", plot_data$condition_nostrain_clean)
plot_data$condition_nostrain_clean = gsub("Pro_PAT_D4", "PAT D4", plot_data$condition_nostrain_clean)


unique(plot_data$sample)

p_ind = ggplot(plot_data, aes(x=condition_nostrain_clean, y=RPKM)) + 
  geom_point(aes(color=strain)) + 
  facet_wrap( ~ gene, ncol=8, scale = "free_y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p_ind


## do the same but now plot the mean value and the standard error per gene and 
## per condition
se <- function(x) sd(x)/sqrt(length(x))

## aggregate by gene, condition_nostrain_clean, and strain
summary_data <- plot_data %>%
  group_by(gene, condition_nostrain_clean, strain) %>%
  summarise(
    mean_RPKM = mean(RPKM),
    se_RPKM = se(RPKM),
    .groups = "drop"
  )

## Select only FC / pval columns we need
pval_df <- dplyr::select(df_out, dplyr::contains("_BPK275_vs_BPK026_fdr")) %>%
  dplyr::mutate(gene = rownames(df_out)) %>%
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "condition_nostrain",
    values_to = "pval"
  ) %>%
  dplyr::mutate(
    condition_nostrain= gsub("_BPK275_vs_BPK026_fdr", "", condition_nostrain),
    signif = dplyr::case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE         ~ ""
    )
  )


pval_df$condition_nostrain_clean = gsub("Pro_LOG_post_PAT", "Post PAT", pval_df$condition_nostrain)
pval_df$condition_nostrain_clean = gsub("Pro_LOG", "No PAT D2", pval_df$condition_nostrain_clean)
pval_df$condition_nostrain_clean = gsub("Pro_PAT_D2", "PAT D2", pval_df$condition_nostrain_clean)
pval_df$condition_nostrain_clean = gsub("Pro_PAT_D4", "PAT D4", pval_df$condition_nostrain_clean)



plot_df <- summary_data %>%
  left_join(pval_df, by = c("gene", "condition_nostrain_clean"))

# asterisk_df <- plot_df %>%
#   group_by(gene, condition_nostrain) %>%
#   summarise(
#     signif = unique(signif[signif != ""]),  # keep only non-empty
#     y_pos = max(mean_RPKM + se_RPKM) + 0.05, # place above error bar
#     .groups = "drop"
#   )
# 
# 
# plot_df$condition_nostrain = factor(plot_df$condition_nostrain, 
#                                       levels= c("Pro_LOG","Pro_PAT_D2",  "Pro_PAT_D4", "Pro_LOG_post_PAT"))

plot_df$condition_nostrain_clean = factor(plot_df$condition_nostrain_clean, 
                                    levels= c("No PAT D2","PAT D2","Post PAT"))


# rpkm_plot = p_ind + ggplot(plot_df, aes(x=condition_nostrain, y=mean_RPKM, color=strain, group=strain)) +
#   geom_point(position=position_dodge(width=0.3), size=2) +
#   geom_errorbar(aes(ymin=mean_RPKM-se_RPKM, ymax=mean_RPKM+se_RPKM),
#                 width=0.2,
#                 position=position_dodge(width=0.3)) +
#   geom_text(data=asterisk_df,
#             aes(x=condition_nostrain, y=y_pos, label=signif),
#             color="black",
#             inherit.aes = FALSE,  # use only asterisk_df columns
#             size=4) +
#   facet_wrap(~ gene, ncol=8, scales="free_y") +
#   scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
#   labs(y="RPKM", x="Condition")
# 
# rpkm_plot




### updated RPKM plot combinding bot the mean + se + individual raw points. 

asterisk_df <- plot_data %>%
  group_by(gene, condition_nostrain_clean) %>%
  summarise(
    y_pos = max(RPKM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    plot_df %>%
      distinct(gene, condition_nostrain_clean, signif),
    by = c("gene", "condition_nostrain_clean")
  ) %>%
  filter(signif != "") %>%
  group_by(gene, condition_nostrain_clean) %>%
  mutate(
    y_pos = y_pos + 0.05 * y_pos   # 5% padding
  ) %>%
  ungroup()


number_fig_columns = 5

rpkm_plot = ggplot(plot_df,
                   aes(x = condition_nostrain_clean,
                       y = mean_RPKM,
                       color = strain,
                       group = strain)) +
  
  ## individual RPKM points (RAW DATA)
  geom_point(
    data = plot_data,
    aes(x = condition_nostrain_clean,
        y = RPKM,
        color = strain),
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width  = 0.3
    ),
    alpha = 0.5,
    size = 1.5,
    inherit.aes = FALSE
  ) +
  
  ## mean points
  geom_point(
    position = position_dodge(width = 0.3),
    size = 2.5
  ) +
  
  ## error bars (SE)
  geom_errorbar(
    aes(ymin = mean_RPKM - se_RPKM,
        ymax = mean_RPKM + se_RPKM),
    width = 0.2,
    position = position_dodge(width = 0.3)
  ) +
  
  ## significance asterisks
  geom_text(
    data = asterisk_df,
    aes(x = condition_nostrain_clean,
        y = y_pos,
        label = signif),
    color = "black",
    size = 4,
    inherit.aes = FALSE
  ) +
  
  facet_wrap(~ gene, ncol = number_fig_columns, scales = "free_y") +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.15))
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "RPKM", x = "Condition")

rpkm_plot

rpkm_plot_file = paste0(out_dir, 'rpkm_plot_file_full.png')

if (number_fig_columns == 8) {
  ggsave(file = rpkm_plot_file,
         plot = rpkm_plot,
         dpi = 300, 
         width = 16,
         height = 9)
}else{
  ## A4 format
  ggsave(file = rpkm_plot_file,
         plot = rpkm_plot,
         dpi = 300, 
         width = 9,
         height = 13)
}




# ---- 4.6 count hypothetical proteins ----

fc_columns = grep("_fc$", colnames(df_out), value=TRUE)
cutoff_fdr = 0.05
cutoff_fc = 1.00
cutoff_fc = 0.75
cutoff_fc = 0.58
res_df = data.frame()

for (fc_column in fc_columns) {
  fdr_column = sub("_fc", "_fdr", fc_column)
  
  ## upregulated genes
  df_out_up = df_out[df_out[,fdr_column] <= cutoff_fdr & df_out[,fc_column] >= cutoff_fc,]
  df_out_up_hypo = length(grep("ypothetical", df_out_up$description))
  
  ## downregulated genes
  df_out_down = df_out[df_out[,fdr_column] <= cutoff_fdr & df_out[,fc_column] <= -cutoff_fc,]
  df_out_down_hypo = length(grep("ypothetical", df_out_down$description))
  
  ## print output
  comparison = sub("_fc", "", fc_column)
  print(paste(comparison, 
              nrow(df_out_up), 
              df_out_up_hypo, 
              nrow(df_out_down), 
              df_out_down_hypo, 
              sep = ","))
  
  res_df <- rbind(
    res_df,
    data.frame(
      comparison = comparison,
      n_up = nrow(df_out_up),
      n_up_hypo = df_out_up_hypo,
      perc_up_hypo = df_out_up_hypo/nrow(df_out_up),
      n_down = nrow(df_out_down),
      n_down_hypo = df_out_down_hypo,
      perc_down_hypo = df_out_down_hypo/nrow(df_out_down),
      stringsAsFactors = FALSE
    )
  )
}

res_df
write.xlsx(res_df, file = paste0(out_dir, 'hypothetical_proteins_counts_fc', cutoff_fc, '_pval', cutoff_fdr, '.xlsx'))

# ---- 4.7 Find those genes uniquely upregulted at PAT exposure ----

## find those genes which are only differentially up (or down) regulated when 
## comparing PAT exposed samples between BPK026 and BPK275, and which are *not*
## differentially expressed when comparing BPK026 and BPK275 in prolog condition
## (either before or after PAT pressure)
head(df_out)
cutoff_fdr = 0.05
cutoff_fc = 1

## for upregulated genes
prolog_up = rownames(df_out[df_out[,"Pro_LOG_BPK275_vs_BPK026_fdr"] <= cutoff_fdr & df_out[,"Pro_LOG_BPK275_vs_BPK026_fc"] >= cutoff_fc,])
postpat_up = rownames(df_out[df_out[,"Pro_LOG_post_PAT_BPK275_vs_BPK026_fdr"] <= cutoff_fdr & df_out[,"Pro_LOG_post_PAT_BPK275_vs_BPK026_fc"] >= cutoff_fc,])
nopat_up = union(prolog_up, postpat_up)
pat_up = rownames(df_out[df_out[,"Pro_PAT_D2_BPK275_vs_BPK026_fdr"] <= cutoff_fdr & df_out[,"Pro_PAT_D2_BPK275_vs_BPK026_fc"] >= cutoff_fc,])

intersect(nopat_up, pat_up)
pat_specific_up = setdiff(pat_up, nopat_up)

df_pat_up = df_out[rownames(df_out) %in% pat_specific_up,]


## run enrichment analysis on KEGG pathways. Run first the code in section 
## 3.8 to prepare everything for the KEGG analysis
ego <- enricher(
  gene = pat_specific_up,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pAdjustMethod = "BH", # BH
  qvalueCutoff = .5,
  pvalueCutoff = 0.05,
  minGSSize = 2
)


enrich_df_pat_specific_up <- as.data.frame(ego)
dim(enrich_df_pat_specific_up)


## for downregulated genes
prolog_down = rownames(df_out[df_out[,"Pro_LOG_BPK275_vs_BPK026_fdr"] <= cutoff_fdr & df_out[,"Pro_LOG_BPK275_vs_BPK026_fc"] <= - cutoff_fc,])
postpat_down = rownames(df_out[df_out[,"Pro_LOG_post_PAT_BPK275_vs_BPK026_fdr"] <= cutoff_fdr & df_out[,"Pro_LOG_post_PAT_BPK275_vs_BPK026_fc"] <= -cutoff_fc,])
nopat_down = union(prolog_down, postpat_down)
pat_down = rownames(df_out[df_out[,"Pro_PAT_D2_BPK275_vs_BPK026_fdr"] <= cutoff_fdr & df_out[,"Pro_PAT_D2_BPK275_vs_BPK026_fc"] <= -cutoff_fc,])

intersect(nopat_down, pat_down)
pat_specific_down = setdiff(pat_down, nopat_down)

df_pat_down = df_out[rownames(df_out) %in% pat_specific_down,]


## run enrichment analysis on KEGG pathways. Run first the code in section 
## 3.8 to prepare everything for the KEGG analysis
ego <- enricher(
  gene = pat_specific_down,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pAdjustMethod = "BH", # BH
  qvalueCutoff = .5,
  pvalueCutoff = 0.05,
  minGSSize = 2
)


enrich_df_pat_specific_down <- as.data.frame(ego)
dim(enrich_df_pat_specific_down)


## create a new workbook
wb <- createWorkbook()

## add empty working sheets
addWorksheet(wb, "PAT_up_genes")
addWorksheet(wb, "PAT_up_enrichment")
addWorksheet(wb, "PAT_down_genes")
addWorksheet(wb, "PAT_down_enrichment")

## write data frames to sheets
writeData(wb, sheet = "PAT_up_genes", df_pat_up)
writeData(wb, sheet = "PAT_up_enrichment", enrich_df_pat_specific_up)
writeData(wb, sheet = "PAT_down_genes", df_pat_down)
writeData(wb, sheet = "PAT_down_enrichment", enrich_df_pat_specific_down)

# Save workbook
saveWorkbook(wb, file = paste0(out_dir, "PAT_specific_genes_and_enrichment.xlsx"), overwrite = TRUE)



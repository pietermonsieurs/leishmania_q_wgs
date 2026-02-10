library(Mfuzz)
library(reshape)
library(pheatmap)
library(openxlsx)
library(ggplot2)
library(Biobase)



#### input data ####



# Variance stabilizing transformation -- start from the bulk_DESeq2.R script
# vsd <- vst(dds, blind = FALSE)   # blind=FALSE uses your design
vsd <- readRDS("/Users/pmonsieurs/programming/leishmania_10X_Q/results/bulk/bwa_htseq/vsd_object_for_mfuzz.R")
expr_matrix <- assay(vsd)        # get the expression matrix

## select the strain for which you want to run the mfuzz clustering. BPK275 
## contains most of the time points. Also, create a vector that contains the 
## different conditions (= time points) that should be used to average the data. 
## For strain BPK275, you can optionally split up D2 samples into D2_early for 
## R1 and R2, and D2_late for R3 and R4

strain = "BPK275"
strain = "BPK026"

expr_matrix_subset = expr_matrix[,grep(strain, colnames(expr_matrix))]
time_points = samples_clean <- sub("_R[0-9]+$", "", colnames(expr_matrix_subset))

## for BPK275 only
time_points[grep("D2",time_points)] = c("BPK275_Pro_PAT_D2_early",
                                        "BPK275_Pro_PAT_D2_early",
                                        "BPK275_Pro_PAT_D2_late",
                                        "BPK275_Pro_PAT_D2_late")

## make factors of them so you can preserve the (chronological) order of the 
## different conditons when calculating the mean expression
time_points <- factor(time_points, levels = unique(time_points))

## average replicates by time point
expr_matrix_mean <- sapply(levels(time_points), function(tp) {
                            idx <- which(time_points == tp)
                            rowMeans(expr_matrix[, idx, drop = FALSE])
                          })

head(expr_matrix_mean)

## Create ExpressionSet, needed as input for the mfuzz, and do stndardization
## of the genes (important for clustering)
eset <- ExpressionSet(assayData = expr_matrix_mean)
eset <- filter.NA(eset, thres=0.01) ## not needed - no NA values present
eset <- filter.std(eset,min.std=0.50)

eset <- standardise(eset)
eset

m <- mestimate(eset)
m

## determination of the number of clusters, only a rough estimate. Best to 
## choose a number of cluster where the min centroid distance goes to a plateau
## so best way is the one where you find the curve in the elbow.
## in the example below, this looks like ~ 8
dmin_results = Dmin(eset, m=m, crange=2:20, repeats=3, visu=TRUE)
dmin_df = as.data.frame(dmin_results)
dmin_df$cluster_number = seq(1,length(dmin_results),1)
colnames(dmin_df)[1] = "min_centroid_distance"

ggplot(data = dmin_df, aes(x=cluster_number, y=min_centroid_distance)) + 
  geom_point(color="blue", size=2) + geom_line(color="darkblue") + 
  theme_bw()


## perform the clustering, in this case with 8 cluster
c <- 5
cl <- mfuzz(eset, c=c, m=m)

## visualize the clustering - with default viz script of mfuzz, but better to 
## reimplement using ggplot?
# mfuzz.plot(eset, cl=cl, mfrow=c(2,3))  # 2 rows × 3 cols

## create a data frame that you can use to plot the same figures as what is
## done using the XQuartz approach
cluster_assign <- data.frame(
  gene = rownames(exprs(eset)),
  cluster = cl$cluster
)

## extract the (standardised) expression information
expr_df <- as.data.frame(exprs(eset))
expr_df$gene <- rownames(expr_df)

## merge the cluster information and the expression information and melt to 
## obtain a dataframe that can be used the make a plot similar as with the
## approach of mfuzz
expr_df <- merge(expr_df, cluster_assign, by="gene")
expr_melt <- melt(expr_df, id.vars = c("gene","cluster"),
                  variable.name = "timepoint", value.name = "expression")
colnames(expr_melt) = c('gene', 'cluster', 'timepoint', 'expression')

ggplot(expr_melt, aes(x = timepoint, y = expression, group = gene)) +
  geom_line(alpha = 0.3) +            # individual gene trends
  stat_summary(aes(group=cluster), 
               fun = mean, geom="line", color="red", linewidth=1.2) +  # cluster mean
  facet_wrap(~ cluster, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


## prepar the export of the full membership and dataframe. The assigned cluster
## ID is the cluster with the highest probability in the membership dataframe. 
## howvever, it can be that this membership is lower than 0.50, simply the highest
## membership value is choosen
mfuzz_df = cbind.data.frame(cl$cluster, cl$membership)
colnames(mfuzz_df)[1] = 'mfuzz_cluster'

## add the gff annotation to the mfuzz_df. Repeat the code from the bulk_DESeq2.R
## to parse the gff file
gff_file = '/Users/pmonsieurs/programming/leishmania_10X_Q/data/refgenome/Leishmania_donovani_16Nov2015beta.gff3'
gff = read.csv(gff_file, sep="\t", comment.char="#", header=FALSE)
unique(gff$V3)
colnames(gff)

## extract the function of a gene
# gff_function = gff[gff$V3 == "protein_coding_gene",]
gff_function = gff
gff_function$description = ""
gff_function$name = ""
#gff_function$geneID = ""

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
gff_sub = gff_function[,c(1,3,4,5,7,10,11,12,13)]
colnames(gff_sub)[1:5] = c('chrom', 'type', 'start', 'stop','strand')
head(gff_sub)


## add the version 1 gene ID
gene_id_file = '/Users/pmonsieurs/programming/leishmania_10X/data/blast/Mapping_Sanger_vs_TriTrypDB.csv'
gene_data = read.table(gene_id_file)
colnames(gene_data) = c('id_v2', 'id_v1')
gene_data$id_v1 <- sub("\\.1$", "", gene_data$id_v1)
gene_data$id_v2 <- sub("\\.1$", "", gene_data$id_v2)

gff_sub$id_v1 = gene_data[match(gff_sub$gene_id, gene_data$id_v2),]$id_v1
mfuzz_df = cbind.data.frame(mfuzz_df, gff_sub[match(rownames(mfuzz_df), gff_sub$gene_id),])

mfuzz_dir = "/Users/pmonsieurs/programming/leishmania_10X_Q/results/bulk/bwa_htseq/mfuzz/"
excel_file = paste0(mfuzz_dir, 'mfuzz.xlsx')
write.xlsx(mfuzz_df, file=excel_file, rowNames=TRUE)


## export gene sets to perform functional enrichment analysis on those gene lists. 
for (cluster in unique(mfuzz_df$mfuzz_cluster)) {
  gene_file = paste0(mfuzz_dir, 'genes_mfuzz_cluster_', cluster, '.txt')
  genes = rownames(mfuzz_df[mfuzz_df$mfuzz_cluster == cluster,])
  write.table(genes, gene_file, quote=FALSE, row.names = FALSE, col.names=FALSE)
}




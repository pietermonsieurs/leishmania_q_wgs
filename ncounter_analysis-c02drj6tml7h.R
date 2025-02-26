## install tools required for doing the analysis of the nCounter 
## independt of the nSolver software
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("NanoStringQCPro")
# BiocManager::install("NanoStringDiff")


# Load the packages
library(NanoStringQCPro)
library(NanoStringDiff)
library(limma)
library(pheatmap)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(openxlsx)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)



#### 1. input data and metadata ####

# ---- ├ 1.1 parameter settings ----
out_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/ncounter/'

## check the meta data. In contrast to what the manual says, the data should
## be stored in a tab-delimited file instead of a comma separated file, so 
## use the .txt metadata and not the .csv version
meta_data_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/'
meta_data_file = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/ncounter_metadata.xlsx'
meta_data = read.xlsx(meta_data_file)
head(meta_data)


# ---- ├ 1.2. nanostring data -----

## specity the rlf file, which is downloaded from the nanostring website
## https://nanostring.com/products/ncounter-assays-panels/immunology/host-response/
rlf_file = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/ncounter/NS_Mm_HostResponse_v1.0.rlf'

rcc_directory <- "/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/ncounter/combined/"
rcc_data <- newRccSet(rccFiles = dir(rcc_directory, full.names = TRUE),
                      extraPdata = meta_data_file,
                      rlf = rlf_file,
                      blankLabel = "blank"
                      )
checkRccSet(rcc_data)

## inspect the resulting data set
str(max.level=2, example_rccSet)
rcc_data@phenoData@data



#### 2. normalisation and QC ####

# ---- ├ 2.1 normaliation ----

rcc_data_norm <- preprocRccSet(rccSet = rcc_data, 
                               normMethod = "housekeeping",
                               bgReference="negatives")


# ---- ├ 2.2 QC function ----

rcc_data_norm_qc <- makeQCReport(rcc_data_norm, 
                                 paste0(out_dir, "nCounter_QC_report.html"),
                                 sampleNameCol = "sample_id")



# ---- ├ 2.3 correlation analysis -----

## correlation between biological replicates
meta_data$group = paste0(meta_data$timepoint, meta_data$infection)
expr_matrix = rcc_data_norm@assayData$normData
colnames(expr_matrix) = unname(sapply(colnames(expr_matrix), function(x) unlist(strsplit(x, split = "\\_"))[3]))
colnames(expr_matrix)[24] = "Standard2"

# filter per "biological condition", e.g. per time point, or timepoint x infection
for (group in unique(meta_data$group)) {
  samples = colnames(expr_matrix)[grep(group, colnames(expr_matrix))]
  print(samples)
}



## calculate the correlation between all samples, and plot the correlations
## in a heatmap. Remove the samples containing the "standard" and optionally
## remove 96hI2 as this seems to be an outlier. Also 24hI2 is present in both
## batches, so one should be renamed (hard coded) to 24hI2.2 (this is also the
## name in the metadata file)
expr_matrix = rcc_data_norm@assayData$normData
colnames(expr_matrix) = unname(sapply(colnames(expr_matrix), function(x) unlist(strsplit(x, split = "\\_"))[3]))
colnames(expr_matrix)[24] = "Standard2"
colnames(expr_matrix)[grep("24hI2", colnames(expr_matrix))[2]] = "24hI2.2"
expr_matrix = expr_matrix[, -grep("Standard", colnames(expr_matrix))]
expr_matrix = expr_matrix[, - which(colnames(expr_matrix) == "96hI2")]
cor_matrix <- cor(expr_matrix, method = "pearson")

annotation_data = meta_data[match(colnames(cor_matrix), meta_data$sample_id),c('timepoint', 'infection', 'batch')]
rownames(annotation_data) = colnames(cor_matrix)

pheatmap(cor_matrix,
         clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         main = "Sample Correlation Heatmap",
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Choose a color palette
         show_rownames = TRUE,  
         show_colnames = TRUE,
         annotation_col = annotation_data
)

grep("Standard", colnames(cor_matrix))


## calculate the R^2 value 

rsq <- function (x, y) cor(x, y) ^ 2
samples = colnames(expr_matrix)
all_plots = list()

## Create an empty grid to store the plots
grid <- expand.grid(i = 1:24, j = 1:24)
count = 0
for (i in 1:(length(samples))) {
  ## first method
  sample1 = samples[i]

  
  for (j in 1:length(samples)) {
    ## second method
    sample2 = samples[j]
    
    # if (i == j) {next}
    
    r2 = rsq(expr_matrix[, sample1], expr_matrix[,sample2])
    
    count = count + 1
    print(paste0("[[", count, "]] ", i, " -- ", j, " = ", r2))
    
    plot_data = as.data.frame(expr_matrix[,c(sample1, sample2)])
    colnames(plot_data) = c('sample1', 'sample2')
    
    p = ggplot(plot_data, aes(x = sample1, y = sample2))    
    p = p + geom_point(alpha=0.3, size=1.5, color="black")
    p = p + geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) 
    p = p + theme_bw()
    # p = p + stat_cor(aes(label = ..rr.label..), color = "red", geom = "label")
    p = p + stat_cor(method = "pearson", label.x = 0, label.y = 5, color="red")
    # p =    
    p
    
    all_plots[[paste0(i, "_", j)]] = p
  }

  
}

p = wrap_plots(all_plots, ncol=24)
correlation_file = paste0(out_dir, 'correlation_samples.png')
ggsave(p, file=correlation_file, width=16, height=9)



# ---- ├  2.4 PCA ----

## do PCA analysis
expr_matrix_scaled <- scale(expr_matrix)
pca_result <- prcomp(t(expr_matrix_scaled))  # Transpose to treat samples as rows
summary(pca_result)  

## add metadata for colouring in PCA plot 
meta_data_ordered <- meta_data[match(rownames(pca_data), meta_data$sample_id), c('timepoint', 'infection', 'batch')]
pca_data <- as.data.frame(pca_result$x)  
pca_data <- cbind.data.frame(pca_data, meta_data_ordered)
pca_data$batch = as.factor(pca_data$batch)

ggplot(pca_data, aes(x = PC1, y = PC2)) +
  # geom_point(aes(color = timepoint), size = 3) +
  # geom_point(aes(color = batch), size = 3) +
  geom_point(aes(color = infection), size = 3) +
  theme_minimal() +
  labs(title = "PCA of Samples", x = "PC1", y = "PC2")


#### 3. differential expression ####

# ---- ├  3.1 using NanoStringDiff ---- 

## obsolete: glm.LRT gets stuck in an endless loop
## do differential expression analysis using NanoStringDiff. For this 
## you need to split up between endogenous probes, positve, negative 
## and housekeeping probbes. 
ncounter_data = rcc_data@assayData[["exprs"]]
endogenous = ncounter_data[grep("Endogenous", rownames(ncounter_data)),]
negative = ncounter_data[grep("Negative", rownames(ncounter_data)),]
positive = ncounter_data[grep("Positive", rownames(ncounter_data)),]
housekeeping = ncounter_data[grep("Housekeeping", rownames(ncounter_data)),]

designs = data.frame(time_point = meta_data$timepoint,
                     infection = meta_data$infection, 
                     group = paste0(meta_data$timepoint, "_", meta_data$infection))
                     
                       
ns_data = createNanoStringSet(endogenous,
                    positive,
                    negative,
                    housekeeping,
                    designs)


pData(ns_data)
head(exprs(ns_data))

## do normalisation on the nanostring data
ns_data = estNormalizationFactors(ns_data)

## create design matrix to perform differential expression analysis
pheno = pData(ns_data)
group = pheno$group
design.full = model.matrix(~0+group)

## do differential expression. This takes forever, somewhere stuck in 
## an endless loop? 
result_12hI =glm.LRT(ns_data, design.full, 
                     Beta=ncol(design.full), 
                     contrast=c(-1,1,0,0,0,0,0,0,0))


# ---- ├ 3.2 edgeR ----

## edgeR on normalised data, but edgeR does not except negative
## values
group <- designs$group # Adjust to your conditions
dge <- DGEList(counts = rcc_data_norm@assayData$normData, group = designs$group)

## Skip TMM normalization for already normalized data. normalised using the 
## NanoStringQCPro package
dge$samples$norm.factors <- 1

design <- model.matrix(~ 0 + group)
dge <- estimateDisp(dge, design)


# ---- ├ 3.3 limma ----

## configure the experimental design
designs = data.frame(time_point = meta_data$timepoint,
                     infection = meta_data$infection, 
                     group = paste0(meta_data$timepoint, "_", meta_data$infection))

g_ = factor(designs$group)
g_ = relevel(g_, ref= "0h_I")
design <- model.matrix(~ 0 + g_)


## fit the linear model
fit <- lmFit(rcc_data_norm@assayData$normData, design)
fit <- eBayes(fit)
fit

## make the contrasts, each time comparing with the control at 0h
contrast_matrix <- makeContrasts(g_4h_I-g_0h_I,
                                g_4h_NI-g_0h_I,
                                g_12h_I-g_0h_I,
                                g_12h_NI-g_0h_I,
                                g_24h_I-g_0h_I,
                                g_24h_NI-g_0h_I,
                                g_96h_I-g_0h_I,
                                g_96h_NI-g_0h_I,
                                levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

fit2$coefficients
fit2$p.value
dim(fit2$coefficients)
dim(fit2$p.value)


## create heatmaps on all genes
pheatmap(fit2$coefficients)

## create heatmap on diff expressed genes
diff_expressed = topTable(fit2, p.value=0.01, lfc=1, n=200) 
rownames(diff_expressed) = gsub("Endogenous_", "", rownames(diff_expressed))
pheatmap(diff_expressed[,1:8],
         fontsize_row = 7, 
         treeheight_col = 0)


## check with output Emma using nSolver. She compared for time point 12h and 
## time point 4h the infected with the non-infected samples
contrast_matrix_emma <- makeContrasts(g_4h_I-g_4h_NI,
                                     g_12h_I-g_12h_NI,
                                     levels = design)
fit_emma <- contrasts.fit(fit, contrast_matrix_emma)
fit_emma <- eBayes(fit_emma)
topTable(fit_emma, coef=2)
topTable(fit_emma, coef=1)



#### 4. GSEA analysis ####

# ---- 4.1 using publicly available databases ----

## Extract gene symbols and convert to EntrezIDS

for (coeff in c(1,3,5,7)) {
  sample_name = colnames(fit2$contrasts)[coeff]
  sample_name = gsub("g_", "", sample_name)
  sample_name = gsub(" ", "", sample_name)
  print(sample_name)
  
  diff_expressed = topTable(fit2, coef = coeff, number = 10000,  p.value=0.05, lfc=1)
  gene_symbols <- sapply(strsplit(rownames(diff_expressed), "_"), `[`, 1)
  entrez_ids <- mapIds(org.Mm.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  print(entrez_ids)
  
  ## GO enrichment analysis
  go_enrich <- enrichGO(gene = entrez_ids,
                        OrgDb = org.Mm.eg.db,
                        ont = "ALL",  ## vary over BP, CC, or MF
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,
                        readable = TRUE,
                        minGSSize = 20)
  
  go_gse <- gseGO(gene = geneList,
                  OrgDb = org.Mm.eg.db,
                  ont = "ALL",  ## vary over BP, CC, or MF
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05)
  dim(go_gse)
  
  ## filter the go_enrich output
  go_enrich_simplified <- simplify(go_enrich, cutoff = 0.5, by = "p.adjust", select_fun = min)
  dim(go_enrich_simplified)
  
  go_gse_simplified <- simplify(go_gse, cutoff = 0.5, by = "p.adjust", select_fun = min)
  dim(go_gse_simplified)
  
  
  
  
  # Visualize GO enrichment results
  go_plot_enriched = dotplot(go_enrich_simplified,
                             showCategory = 25,
                             font.size = 7,
  )
  
  
  go_plot_gse = dotplot(go_gse_simplified,
                        showCategory = 25,
                        font.size = 10,
  )
  
  go_plot_file_enriched = paste0(out_dir, 'gsea/', sample_name, '_enricheGO.png')
  go_plot_file_gse = paste0(out_dir, 'gsea/', sample_name, '_gsego.png')
  
  ggsave(go_plot_enriched, file = go_plot_file_enriched)
  ggsave(go_plot_gse, file = go_plot_file_gse)

}


  
  
  # ---- 4.2 using categories of Nanostring using hypergeometric distribution ----
  

## To Do ##
## correlation between normalised values of the biological replicates
## enhanced volcanoplot per condition






## load normalised data into 



#### example data & code snippets ####
exampleDataDir <- system.file("extdata", package="NanoStringQCPro")
rccDir <- file.path(exampleDataDir, "RCC")
example_rccSet <- newRccSet(
  rccFiles = dir(rccDir, full.names=TRUE)
  ,rlf = file.path(exampleDataDir, "RLF", "NQCP_example.rlf")
  ,cdrDesignData = file.path(exampleDataDir, "CDR", "CDR-DesignData.csv")
  ,extraPdata = file.path(exampleDataDir, "extraPdata", "SampleType.txt")
  ,blankLabel = "blank"
)




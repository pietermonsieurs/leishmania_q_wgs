## install tools required for doing the analysis of the nCounter 
## independt of the nSolver software
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("NanoStringQCPro")
# BiocManager::install("NanoStringDiff")


# Load the packages
library(NanoStringQCPro)
library(NanoStringDiff)
library(edgeR)
library(pheatmap)




#### 1. input data and metadata ####

# ---- 1.0 parameter settings ----
out_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/ncounter/'

# ---- 1.1. nanostring data -----

## check the meta data. In contrast to what the manual says, the data should
## be stored in a tab-delimited file instead of a comma separated file, so 
## use the .txt metadata and not the .csv version
meta_data_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/'
meta_data_file = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/ncounter_metadata.txt'

meta_data = read.table(meta_data_file, header=TRUE)

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



rcc_data_norm <- preprocRccSet(rccSet = rcc_data, 
                               normMethod = "housekeeping",
                               bgReference="negatives")

rcc_data_norm_qc <- makeQCReport(rcc_data_norm, 
                                 paste0(out_dir, "example_QC_report"))



## do differential expression analysis. For this you need to split up between
## endogenous probes, positve, negative and housekeeping probbes
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

## do differential expression 
result_12hI =glm.LRT(ns_data, design.full, Beta=ncol(design.full), contrast=c(-1,1,0,0,0,0,0,0,0))



## edgeR 
group <- designs$group # Adjust to your conditions
dge <- DGEList(counts = rcc_data_norm@assayData$normData, group = designs$group)


## Skip TMM normalization for already normalized data. normalised using the 
## NanoStringQCPro package
dge$samples$norm.factors <- 1

design <- model.matrix(~ 0 + group)
dge <- estimateDisp(dge, design)

group = factor(designs$group)
group = relevel(group, ref= "0h_I")
design <- model.matrix(~ 0+ group)


# Fit the linear model
fit <- lmFit(rcc_data_norm@assayData$normData, design)
fit <- eBayes(fit)
fit

contrast_matrix <- makeContrasts(group4h_I-group0h_I,
              group4h_NI-group0h_I,
              group12h_I-group0h_I,
              group12h_NI-group0h_I,
              group24h_I-group0h_I,
              group24h_NI-group0h_I,
              group96h_I-group0h_I,
              group96h_NI-group0h_I,
              levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

fit2$coefficients
dim(fit2$coefficients)

pheatmap(fit2$coefficients)

diff_expressed = topTable(fit2, p.value=0.01, lfc=1, n=200) 


diff_expressed = topTable(fit2, p.value=0.01, lfc=1, n=200) 
pheatmap(diff_expressed[,1:8],
         fontsize_row = 6)


## check with output Emma using nSolver. She compared for time point 12h and 
## time point 4h the infected with the non-infected samples
contrast_matrix_emma <- makeContrasts(group4h_I-group4h_NI,
                                 group12h_I-group12h_NI,
                                 levels = design)
fit_emma <- contrasts.fit(fit, contrast_matrix_emma)
fit_emma <- eBayes(fit_emma)
topTable(fit_emma, coef=2)

## correlation between normalised values of the biological replicates? 





## load normalised data into 



#### example data ####
exampleDataDir <- system.file("extdata", package="NanoStringQCPro")
rccDir <- file.path(exampleDataDir, "RCC")
example_rccSet <- newRccSet(
  rccFiles = dir(rccDir, full.names=TRUE)
  ,rlf = file.path(exampleDataDir, "RLF", "NQCP_example.rlf")
  ,cdrDesignData = file.path(exampleDataDir, "CDR", "CDR-DesignData.csv")
  ,extraPdata = file.path(exampleDataDir, "extraPdata", "SampleType.txt")
  ,blankLabel = "blank"
)

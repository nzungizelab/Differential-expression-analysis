
# set working directory
setwd("D:/ML/K-mean/jako")

# path for saving ouptut
outpath <- "D:/ML/K-mean/jako/"  

#loading data (counts matrix)
counts <-read.delim("joka_normalization.txt", row.names = 1, sep="\t", header = TRUE) # counts can be read into a data.frame

# view of data loaded
head(counts)
dim(counts)
colnames(counts)
data(counts)

# load genotypes data and conditions 
expdesign1 <- read.table("pheno.txt", header=T)
head(expdesign1, n = 15)
dim(expdesign1)

# check the groups we have in our experiment design
groups <- paste(expdesign1$strain2, expdesign1$stage2, sep=".")
groups <- factor(groups) #group	vector containing the experimental group/condition for each sample(library)
table(groups) 

# load libraries
library("Rcpp")
library("NMF")
library("DESeq2")
library("edgeR")
library("limma")
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("EnhancedVolcano")
library("pheatmap")

# create DGEList object using reads counts 
class(counts) #checking if is data frame
counts <- as.matrix(counts) #convert df to matrix, because DGEList needs the counts as a matrix 

y <- DGEList(counts, samples=expdesign1, group=expdesign1$strain2)
head(y)

# Filtering data
#using edgeR the raw counts are converted to CPM and log-CPM values using the cpm function.
#By default in edgeR, a gene should have CPM > 0.5 in at least 3 samples in opposite condition, the gene is removed

dim(y) #Checking before filtering step
#[1] 4681   22

keep <- rowSums(cpm(y) > 0.5) >= 2 #Identify and only keep genes with cpm value greater than 0.5 in at least 2 samples.

table(keep)

y_filtered <- y[keep, , keep.lib.sizes = FALSE]

head(y_filtered)

dim(y_filtered) #Checking after filtering step

#Normalization
#Among various normalisation methods (TMM, BH, RLE, and upperquartile UQ) we will use  TMM (most useful normalization in RNA-Seq analysis)
#These normalization methods are performed using either DESeq (RLE) or edgeR (TC, RPKM, UQ, TMM).
#TMM(trimmed mean of M-values: estimating relative RNA production levels from RNA-seq data

y_filtered <- calcNormFactors(y_filtered, method="TMM")  #Calculate Normalization Factors to Align Columns of a Count Matrix 
head(y_filtered$samples, n=15)

#creating matrix
#design matrix helps to compare groups of your data
design <- model.matrix(~0+groups)

colnames(design) <- levels(groups)
dim(design)
design

y_filtered <- estimateDisp(y_filtered, design, robust=TRUE)
y_filtered$common.dispersion
#[1] 9.765625e-05

#edgeR uses the Cox-Reid profile-adjusted likelihood method in estimating dispersion
y_filtered <- estimateCommonDisp(y_filtered,verbose=TRUE)
#Disp = 1e-04 , BCV = 0.01 

#Estimated values of the GLM coefficients for each gene.
fit <- glmQLFit(y_filtered, design, robust=TRUE) 

head(fit$coefficients)
summary(fit$df.prior)

#After filtering and normalization, check how samples are relate to each other 
plotMDS(y_filtered, method="bcv", col=as.numeric(y_filtered$samples$groups))
legend("topleft", as.character(unique(y_filtered$samples$groups)), col=1:3, pch=20)



# comparing two groups to finds out upregulated and downregulated genes 
# comparing 2 groups we can contrast() using the method exactTest(DGEList) to obtain DEG. 
#more than 2 groups use a linear model. For edgeR uses a generalized linear model ( GLM )

labSexvsAsex <- makeContrasts(Lab.Asexual-Lab.Sexual, levels=design)  

clinicalSexvsAsex <- makeContrasts(Clinical.Asexual-Clinical.Sexual, levels=design)

sex <- makeContrasts(Clinical.Sexual-Lab.Sexual, levels=design)

asexual <- makeContrasts(Clinical.Asexual-Lab.Asexual, levels=design)



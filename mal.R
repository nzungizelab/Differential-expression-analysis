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


#load gene annotation
annot <-read.delim("gene_ane_symbols1.txt", header = TRUE) # counts can be read into a data.frame

# view of data loaded
head(annot)
dim(annot)
colnames(annot)

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

##Data exploration

# using PCA
#After filtering and normalization, check how samples are relate to each other 
png(paste0(outpath,"PCA1.png"))
plotMDS(y_filtered, 
        col=c(rep("black",9), rep("red",7), rep("green",3),rep("yellow",3)), 
        labels= c(rep("Clinical.Asexual",9), rep("Clinical.Sexual",7), rep("Lab.Asexual",3),rep("Lab.Sexual",3)))

dev.off()

#chek sample data without group label
png(paste0(outpath,"PCA2.png"))
plotMDS(y_filtered, 
        col=c(rep("black",9), rep("red",7), rep("green",3),rep("yellow",3)))

dev.off()

#using heatmap
#Plotting correlation of samples with normalized counts
png(paste0(outpath,"correlation.png"))
pheatmap(cor(log2(cpm(y_filtered)+1), cluster_cols = F,annotation_row = y_filtered), main="Correlation")
dev.off()

# Comparing two groups to finds out upregulated and downregulated genes 
# by comparing 2 groups we can contrast() using the method exactTest(DGEList) to obtain DEG. 
# more than 2 groups use a linear model. For edgeR uses a generalized linear model ( GLM )

labSexvsAsex <- makeContrasts(Lab.Asexual-Lab.Sexual, levels=design)  

clinicalSexvsAsex <- makeContrasts(Clinical.Asexual-Clinical.Sexual, levels=design)

sex <- makeContrasts(Clinical.Sexual-Lab.Sexual, levels=design)

asexual <- makeContrasts(Clinical.Asexual-Lab.Asexual, levels=design)


## Differential Expressed in each group (using the QL F-test).

# 1. Lab strain for Sexual versus Asexual stage
Tab1 <- glmQLFTest(fit, contrast=labSexvsAsex, coef = 2)

topTags(Tab1, n= 5) # table of top 5 DEG

Tab1_top.gene <- topTags(Tab1)  # table of top DEG
Tab1_top.gene

Tab1_de <- decideTestsDGE(Tab1, p.value=0.05) # pvalue using Wald test method
#pvalue=0.01(Down 0, Up 6), pvalue=0.05(Down 7, Up 21)

summary(Tab1_de) 
# Down o gene and Up 0 gene

# adjust p-values and assign the result to our table
#we consider a fraction of 10% false positives acceptable, therefore all genes with an adjusted p value below 10% = 0.1 as significant. 
Tab1$table$padj <- p.adjust(Tab1$table$PValue, method="BH")

sum(Tab1$table$padj < 0.05)
# 28 gene

head(Tab1$table$padj)

topTags(Tab1)

# How many genes look significant?
sum(Tab1$table$PValue < 0.05)
# 397 genes

sum(Tab1$table$padj < 0.05)
# padj < 0.05 = 28

# How many genes show 2-fold enrichment?
sum(Tab1$table$PValue < 0.05 & Tab1$table$logFC > 1)
#99 genes

sum(Tab1$table$padj < 0.05 & Tab1$table$logFC > 1)
# padj < 0.05 = 21


# 2. Clinical strain for Sexual versus Asexual stage
Tab2 <- glmQLFTest(fit, contrast=clinicalSexvsAsex, coef = 2)

topTags(Tab2, n= 5) # table of top 5 DEG

Tab2_top.gene <- topTags(Tab2)  # table of top DEG
Tab2_top.gene

Tab2_de <- decideTestsDGE(Tab2, p.value=0.05) # pvalue using Wald test method
#pvalue=0.01(Down 0, Up 6), pvalue=0.05(Down 7, Up 21)

summary(Tab2_de) 
# Down 195 genes and Up 32 genes

# adjust p-values and assign the result to our table
#we consider a fraction of 10% false positives acceptable, therefore all genes with an adjusted p value below 10% = 0.1 as significant. 
Tab2$table$padj <- p.adjust(Tab2$table$PValue, method="BH")

sum(Tab2$table$padj < 0.05)
# 227 genes

head(Tab2$table$padj)

topTags(Tab2)

# How many genes look significant?
sum(Tab2$table$PValue < 0.05)
# 542 genes

sum(Tab2$table$padj < 0.05)
# padj < 0.05 = 227

# How many genes show 2-fold enrichment?
sum(Tab2$table$PValue < 0.05 & Tab2$table$logFC > 1)
# 8 genes

sum(Tab2$table$padj < 0.05 & Tab2$table$logFC > 1)
# padj < 0.05 = 8


# 3. At Sexual stage between clinical and lab strain
Tab3 <- glmQLFTest(fit, contrast=sex, coef = 2)

topTags(Tab3, n= 5) # table of top 5 DEG

Tab3_top.gene <- topTags(Tab3)  # table of top DEG
Tab3_top.gene

Tab3_de <- decideTestsDGE(Tab3, p.value=0.05) # pvalue using Wald test method
#pvalue=0.01(Down 0, Up 6), pvalue=0.05(Down 7, Up 21)

summary(Tab3_de) 
# Down 0 genes and Up 0 genes

# adjust p-values and assign the result to our table
#we consider a fraction of 10% false positives acceptable, therefore all genes with an adjusted p value below 10% = 0.1 as significant. 
Tab3$table$padj <- p.adjust(Tab3$table$PValue, method="BH")

sum(Tab3$table$padj < 0.05)
# 0 genes

head(Tab3$table$padj)

topTags(Tab3)

# How many genes look significant?
sum(Tab3$table$PValue < 0.05)
# 30 genes

sum(Tab3$table$padj < 0.05)
# padj < 0.05 = 0 gene

# How many genes show 2-fold enrichment?
sum(Tab3$table$PValue < 0.05 & Tab3$table$logFC > 1)
# 11 genes

sum(Tab3$table$padj < 0.05 & Tab3$table$logFC > 1)
# padj < 0.05 = 0 gene


# 4. At Asexual stage between clinical and lab strain
Tab4 <- glmQLFTest(fit, contrast=asexual, coef = 2)

topTags(Tab4, n= 5) # table of top 5 DEG

Tab4_top.gene <- topTags(Tab4)  # table of top DEG
Tab4_top.gene

Tab4_de <- decideTestsDGE(Tab4, p.value=0.05) # pvalue using Wald test method
#pvalue=0.01(Down 0, Up 6), pvalue=0.05(Down 7, Up 21)

summary(Tab4_de) 
# Down 2 genes and Up 0 genes

# adjust p-values and assign the result to our table
#we consider a fraction of 10% false positives acceptable, therefore all genes with an adjusted p value below 10% = 0.1 as significant. 
Tab4$table$padj <- p.adjust(Tab4$table$PValue, method="BH")

sum(Tab4$table$padj < 0.05)
# 2 genes

head(Tab4$table$padj)

topTags(Tab4)

# How many genes look significant?
sum(Tab4$table$PValue < 0.05)
# 50 genes

sum(Tab4$table$padj < 0.05)
# padj < 0.05 =2 gene

# How many genes show 2-fold enrichment?
sum(Tab4$table$PValue < 0.05 & Tab4$table$logFC > 1)
# 3 genes

sum(Tab4$table$padj < 0.05 & Tab4$table$logFC > 1)
# padj < 0.05 = 0 gene


#Volcano plot for Tab1 

#change into data frame
Tab1

str(Tab1) #check data structure

Tab1 <- as.data.frame(Tab1) #save data as df

#Tab1$Gene_Name

head(Tab1)

dim(Tab1)
#[1] 4681    5

#set name of first column header as "GeneID" 
library(dplyr)
library(tibble)

Tab1 <- Tab1 %>% rownames_to_column(var = "GeneID") #Assign a variable to first column

head(Tab1)

#merge two tables (counts and gene annotation table keytype:GeneId)
head(Tab1) #table 1
dim(Tab1)
head(annot) #table 2
dim(annot)

Tab1 <- merge(Tab1,annot, by="GeneID")
head(Tab1, n=10)
dim(Tab1)

#volvano plot
#Volcano plots represent a useful way to visualize the results of DE analyses.
png(paste0(outpath,"volcano_Tab1.png"))
EnhancedVolcano(Tab1,
                #lab = rownames(Tab1),
                lab = Tab1$GeneName,
                x = "logFC",
                y = "padj",
                #y = "PValue",
                #xlim=c(-8, 8),
                ylim=c(0, -log10(10e-6)),
                xlab = bquote(~Log[2]~ 'fold change'),
                #ylab = bquote(~Log[10]~ italic('PValue')),
                ylab = bquote(-~Log[10]~italic('padj')),
                title = 'Lab strain Sexual vs Asexual',
                #selectLab = c('LRP2','TNMD','AFF3','CCDC124'), # lab a set of genes pass FCcutoff and pCutoff thresholds
                
                #labFace = 'bold', #draw labels in bold
                #boxedLabels = TRUE, #draw labels in boxes
                pCutoff = 50e-3, #Cut-off for statistical significance. (horizontal line drawn at -log10(pCutoff),by default = 0.05) also can be (1.3,  10e-6, 0.0001, 10e-2)
                #pLabellingCutoff = pCutoff, #Labelling cut-off for statistical significance. by default = pCutoff
                FCcutoff = 0.5, #Cut-off for absolute log2 fold-change. Vertical lines drawn at the negative/positive values of FCCutoff. by deafault=2.0 also can be (1.3, 1.5,  1.7)
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                #pointSize = 2,
                labSize = 6, #gene label size
                colCustom = NULL,
                #colAlpha = 1,
                #legendLabSize = 15,
                legendPosition = 'bottom',
                legendIconSize = 5.0,
                #drawConnectors = TRUE,
                #widthConnectors = 0.4,
                #cutoffLineType = "longdash",
                #cutoffLineCol = "black",
                #cutoffLineWidth = 0.4,
                #gridlines.major = TRUE,
                #gridlines.minor = TRUE
                #colConnectors = "grey50",
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'partial',
                #borderWidth = 0.5,
                #borderColour = 'black'
                #shape = c(1, 4, 23, 25), #each type of gene expressed with different shape
                #drawConnectors = TRUE,
                #colConnectors = 
                widthConnectors = 0.75,
                #typeConnectors =  # can be "open" or "closed",
                #endsConnectors =   # can be ("last", "first", "both")
                #lengthConnectors = (default = unit(0.01, 'npc'))
                subtitle = NULL,
                #subtitle = "XXXX",
                #caption = paste0("total = ", nrow(toptable), " variables"),
                caption = NULL,
                
)
dev.off()


#Volcano plot for Tab2 

#change into data frame
Tab2

str(Tab2)

Tab2 <- as.data.frame(Tab2)

head(Tab2)

dim(Tab2)

#name column header for first column as "GeneID"
Tab2 <- Tab2 %>% rownames_to_column(var = "GeneID") #Assign a variable to first column

head(Tab2)

#merge two tables (counts and gene annotation table keytype:Geneid)
head(Tab2) #table 1
dim(Tab2)
head(annot) #table 2
dim(annot)

Tab2 <- merge(Tab2,annot, by="GeneID")
head(Tab2, n=10)
dim(Tab2)


#volvano plot
#Volcano plots represent a useful way to visualize the results of differential expression analyses.

png(paste0(outpath,"volcano_Tab2.png"))
EnhancedVolcano(Tab2,
                #lab = rownames(Tab2),
                lab = Tab2$GeneName,
                x = "logFC",
                y = "padj",
                #y = "PValue",
                #xlim=c(-8, 8),
                ylim=c(0, -log10(10e-6)),
                xlab = bquote(~Log[2]~ 'fold change'),
                #ylab = bquote(~Log[10]~ italic('PValue')),
                ylab = bquote(-~Log[10]~italic('padj')),
                title = 'Clinical strain Sexual vs Asexual',
                #selectLab = c('LRP2','TNMD','AFF3','CCDC124'), # lab a set of genes pass FCcutoff and pCutoff thresholds
                
                #labFace = 'bold', #draw labels in bold
                #boxedLabels = TRUE, #draw labels in boxes
                pCutoff = 50e-3, #Cut-off for statistical significance. (horizontal line drawn at -log10(pCutoff),by default = 0.05) also can be (1.3,  10e-6, 0.0001, 10e-2)
                #pLabellingCutoff = pCutoff, #Labelling cut-off for statistical significance. by default = pCutoff
                FCcutoff = 0.5, #Cut-off for absolute log2 fold-change. Vertical lines drawn at the negative/positive values of FCCutoff. by deafault=2.0 also can be (1.3, 1.5,  1.7)
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                #pointSize = 2,
                labSize = 6, #gene label size
                colCustom = NULL,
                #colAlpha = 1,
                #legendLabSize = 15,
                legendPosition = 'bottom',
                legendIconSize = 5.0,
                #drawConnectors = TRUE,
                #widthConnectors = 0.4,
                #cutoffLineType = "longdash",
                #cutoffLineCol = "black",
                #cutoffLineWidth = 0.4,
                #gridlines.major = TRUE,
                #gridlines.minor = TRUE
                #colConnectors = "grey50",
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'partial',
                #borderWidth = 0.5,
                #borderColour = 'black'
                #shape = c(1, 4, 23, 25), #each type of gene expressed with different shape
                #drawConnectors = TRUE,
                #colConnectors = 
                widthConnectors = 0.75,
                #typeConnectors =  # can be "open" or "closed",
                #endsConnectors =   # can be ("last", "first", "both")
                #lengthConnectors = (default = unit(0.01, 'npc'))
                subtitle = NULL,
                #subtitle = "XXXX",
                #caption = paste0("total = ", nrow(toptable), " variables"),
                caption = NULL,
                
)
dev.off()

#Volcano plot for Tab3 

#change into data frame
Tab3

str(Tab3)

Tab3 <- as.data.frame(Tab3)

head(Tab3)

dim(Tab3)


#name column header for first column as "GeneID"
Tab3 <- Tab3 %>% rownames_to_column(var = "GeneID") #Assign a variable to first column

head(Tab3)

#merge two tables (counts and gene annotation table keytype:Geneid)
head(Tab3) #table 1
dim(Tab3)
head(annot) #table 2
dim(annot)

Tab3 <- merge(Tab3,annot, by="GeneID")
head(Tab3, n=10)
dim(Tab3)

#volvano plot
#Volcano plots represent a useful way to visualize the results of differential expression analyses.

png(paste0(outpath,"volcano_Tab3.png"))
EnhancedVolcano(Tab3,
                #lab = rownames(Tab2),
                lab = Tab3$GeneName,
                x = "logFC",
                y = "padj",
                #y = "PValue",
                #xlim=c(-8, 8),
                ylim=c(0, -log10(10e-6)),
                xlab = bquote(~Log[2]~ 'fold change'),
                #ylab = bquote(~Log[10]~ italic('PValue')),
                ylab = bquote(-~Log[10]~italic('padj')),
                title = 'Sexual stage between  Lab and Clinical',
                #selectLab = c('LRP2','TNMD','AFF3','CCDC124'), # lab a set of genes pass FCcutoff and pCutoff thresholds
                
                #labFace = 'bold', #draw labels in bold
                #boxedLabels = TRUE, #draw labels in boxes
                pCutoff = 50e-3, #Cut-off for statistical significance. (horizontal line drawn at -log10(pCutoff),by default = 0.05) also can be (1.3,  10e-6, 0.0001, 10e-2)
                #pLabellingCutoff = pCutoff, #Labelling cut-off for statistical significance. by default = pCutoff
                FCcutoff = 0.5, #Cut-off for absolute log2 fold-change. Vertical lines drawn at the negative/positive values of FCCutoff. by deafault=2.0 also can be (1.3, 1.5,  1.7)
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                #pointSize = 2,
                labSize = 6, #gene label size
                colCustom = NULL,
                #colAlpha = 1,
                #legendLabSize = 15,
                legendPosition = 'bottom',
                legendIconSize = 5.0,
                #drawConnectors = TRUE,
                #widthConnectors = 0.4,
                #cutoffLineType = "longdash",
                #cutoffLineCol = "black",
                #cutoffLineWidth = 0.4,
                #gridlines.major = TRUE,
                #gridlines.minor = TRUE
                #colConnectors = "grey50",
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'partial',
                #borderWidth = 0.5,
                #borderColour = 'black'
                #shape = c(1, 4, 23, 25), #each type of gene expressed with different shape
                #drawConnectors = TRUE,
                #colConnectors = 
                widthConnectors = 0.75,
                #typeConnectors =  # can be "open" or "closed",
                #endsConnectors =   # can be ("last", "first", "both")
                #lengthConnectors = (default = unit(0.01, 'npc'))
                subtitle = NULL,
                #subtitle = "XXXX",
                #caption = paste0("total = ", nrow(toptable), " variables"),
                caption = NULL,
                
)
dev.off()


#Volcano plot for Tab4 

#change into data frame
Tab4

str(Tab4)

Tab4 <- as.data.frame(Tab4)

head(Tab4)
dim(Tab4)

#name column header for first column as "GeneID"
Tab4 <- Tab4 %>% rownames_to_column(var = "GeneID") #Assign a variable to first column
head(Tab4)

#merge two tables (counts and gene annotation table key type:GeneID)
head(Tab4) #table 1
dim(Tab4)
head(annot) #table 2
dim(annot)

Tab4 <- merge(Tab4,annot, by="GeneID")
head(Tab4, n=10)
dim(Tab4)

#volvano plot
#Volcano plots represent a useful way to visualize the results of differential expression analyses.

png(paste0(outpath,"volcano_Tab4.png"))
EnhancedVolcano(Tab4,
                #lab = rownames(Tab4),
                lab = Tab4$GeneName,
                x = "logFC",
                y = "padj",
                #y = "PValue",
                #xlim=c(-8, 8),
                ylim=c(0, -log10(10e-6)),
                xlab = bquote(~Log[2]~ 'fold change'),
                #ylab = bquote(~Log[10]~ italic('PValue')),
                ylab = bquote(-~Log[10]~italic('padj')),
                title = 'Asexual stage between  Lab and Clinical',
                #selectLab = c('LRP2','TNMD','AFF3','CCDC124'), # lab a set of genes pass FCcutoff and pCutoff thresholds
                
                #labFace = 'bold', #draw labels in bold
                #boxedLabels = TRUE, #draw labels in boxes
                pCutoff = 50e-3, #Cut-off for statistical significance. (horizontal line drawn at -log10(pCutoff),by default = 0.05) also can be (1.3,  10e-6, 0.0001, 10e-2)
                #pLabellingCutoff = pCutoff, #Labelling cut-off for statistical significance. by default = pCutoff
                FCcutoff = 0.5, #Cut-off for absolute log2 fold-change. Vertical lines drawn at the negative/positive values of FCCutoff. by deafault=2.0 also can be (1.3, 1.5,  1.7)
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                #pointSize = 2,
                labSize = 6, #gene label size
                colCustom = NULL,
                #colAlpha = 1,
                #legendLabSize = 15,
                legendPosition = 'bottom',
                legendIconSize = 5.0,
                #drawConnectors = TRUE,
                #widthConnectors = 0.4,
                #cutoffLineType = "longdash",
                #cutoffLineCol = "black",
                #cutoffLineWidth = 0.4,
                #gridlines.major = TRUE,
                #gridlines.minor = TRUE
                #colConnectors = "grey50",
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'partial',
                #borderWidth = 0.5,
                #borderColour = 'black'
                #shape = c(1, 4, 23, 25), #each type of gene expressed with different shape
                #drawConnectors = TRUE,
                #colConnectors = 
                widthConnectors = 0.75,
                #typeConnectors =  # can be "open" or "closed",
                #endsConnectors =   # can be ("last", "first", "both")
                #lengthConnectors = (default = unit(0.01, 'npc'))
                subtitle = NULL,
                #subtitle = "XXXX",
                #caption = paste0("total = ", nrow(toptable), " variables"),
                caption = NULL,
                
)
dev.off()

# Venn digram (Tab1 and Tab2)
library("grid")
library("VennDiagram")
library("gplots")

#set value for level of DEG
pval_threshold <- 0.05 
logfc_threshold <- 0.5

#set gene for Tab1 (Lab strain)
Tab1 <- glmQLFTest(fit, contrast=labSexvsAsex, coef = 2)

topTags(Tab1, n= 5) # table of top 5 DEG

# assign adjust p-values to our table
#we consider a fraction of 5% false positives acceptable, therefore all genes with an adjusted p value below 5% = 0.05 as significant. 
Tab1$table$padj <- p.adjust(Tab1$table$PValue, method="BH")

sum(Tab1$table$padj < 0.05)
# 28 gene

head(Tab1$table$padj)

topTags(Tab1)

# Filter on adjusted p-value and get the rownames
Tab1 <- as.data.frame(topTags(Tab1, n = nrow(Tab1)))

#identify the DEG according to the threshold we set
Tab1_Top.deg <- row.names(Tab1[Tab1$FDR <= pval_threshold,])


#set data Tab2 (Clinical strain)
Tab2 <- glmQLFTest(fit, contrast=clinicalSexvsAsex, coef = 2)

topTags(Tab2, n= 5) # table of top 5 DEG

# add adjust p-values 
Tab2$table$padj <- p.adjust(Tab2$table$PValue, method="BH")

sum(Tab2$table$padj < 0.05)
# 227 genes

head(Tab2$table$padj)

topTags(Tab2)

# Filter on adjusted p-value and get the rownames
Tab2 <- as.data.frame(topTags(Tab2, n = nrow(Tab2)))

Tab2_Top.deg <- row.names(Tab2[Tab2$FDR <= pval_threshold,])

# Calculate the intersection of the two sets
deg.intersect = length(intersect(Tab1_Top.deg, Tab2_Top.deg))
deg.venn <- list('intersect' = deg.intersect,
                 'lab' = length(Tab1_Top.deg),
                 'clinical' = length(Tab2_Top.deg))



png(paste0(outpath,"venn_diagramm2.png"))
venn.plot <- draw.pairwise.venn(deg.venn$lab, deg.venn$clinical, deg.venn$intersect,
                                category = c("Lab strain", "Clinical strain"), #label for each circle
                                scaled = F,
                                #filename = '7.venn_diagramm1.png',
                                lwd = rep(2, 2),
                                lty = rep("solid", 2),
                                col = "transparent",
                                fill = c("pink", "green"), #add color to each circles
                                #fill = mypalete3,
                                fontface = "bold",
                                alpha = rep(0.5, 2),
                                cex = rep(1, 3), 
                                cat.pos = c(0, 0),
                                cat.dist = rep(0.025, 2),
                                na = "remove"
)
dev.off()

##Save the DEG into .csv table

#load table with DEG in Lab strain
Tab1_Top.deg <- as.data.frame(Tab1_Top.deg)
head(Tab1_Top.deg)
dim(Tab1_Top.deg)

#rename the header as "GeneID"
Tab1_Top.deg = rename(Tab1_Top.deg, GeneID = Tab1_Top.deg)
head(Tab1_Top.deg)

#creating a table with padj, FDR, LogFc, ...
Tab1 <- glmQLFTest(fit, contrast=labSexvsAsex, coef = 2)

topTags(Tab1, n= 5) # table of top 5 DEG

# add adjust p-values to our table
Tab1$table$padj <- p.adjust(Tab1$table$PValue, method="BH")

sum(Tab1$table$padj < 0.05)
# 28 gene

head(Tab1$table$padj)

topTags(Tab1)

# Filter top DEG based on adjusted p-value and get the rownames
Tab1 <- as.data.frame(topTags(Tab1, n = nrow(Tab1)))
head(Tab1)

#rename the header (first column)
Tab1<- Tab1 %>% rownames_to_column(var = "GeneID") #Assign a variable to first column
head(Tab1)
dim(Tab1)

#merge the two tables (to identify DEG from lab and their padj)
Tab1_Top.deg <- merge(Tab1,Tab1_Top.deg, by="GeneID")
head(Tab1_Top.deg, n=10)
dim(Tab1_Top.deg)


#merge the two tables (DEG gene from Lab strain with annoted table)
Tab1_Top.deg <- merge(Tab1_Top.deg,annot, by="GeneID")
head(Tab1_Top.deg, n=10)
dim(Tab1_Top.deg)

#save the DEG in table as csv format
write.csv(Tab1_Top.deg,"Tab1_Top.deg.csv",row.names=FALSE,quote=FALSE)


#load table with DEG in clinical strain
Tab2_Top.deg <- as.data.frame(Tab2_Top.deg)
head(Tab2_Top.deg)
dim(Tab2_Top.deg)

#rename the header as "GeneID"
Tab2_Top.deg = rename(Tab2_Top.deg, GeneID = Tab2_Top.deg)
head(Tab2_Top.deg)

#creating a table with padj, FDR, LogFc, ...
Tab2 <- glmQLFTest(fit, contrast=clinicalSexvsAsex, coef = 2)

topTags(Tab2, n= 5) # table of top 5 DEG

# adjust p-values and assign the result to our table
Tab2$table$padj <- p.adjust(Tab2$table$PValue, method="BH")

sum(Tab2$table$padj < 0.05)
# 28 gene

head(Tab2$table$padj)

topTags(Tab2)

# Filter on adjusted p-value and get the rownames
Tab2 <- as.data.frame(topTags(Tab2, n = nrow(Tab2)))
head(Tab2)

#rename the header
Tab2<- Tab2 %>% rownames_to_column(var = "GeneID") #Assign a variable to first column
head(Tab2)
dim(Tab2)

#merge the two tables for LogFC, FDR, padj,..
Tab2_Top.deg <- merge(Tab2,Tab2_Top.deg, by="GeneID")
head(Tab2_Top.deg, n=10)
dim(Tab2_Top.deg)

#merge the two tables for LogFC, FDR, padj,.. with annoted table
Tab2_Top.deg <- merge(Tab2_Top.deg,annot, by="GeneID")
head(Tab2_Top.deg, n=10)
dim(Tab2_Top.deg)

#save the table in csv format
write.csv(Tab2_Top.deg,"Tab2_Top.deg.csv",row.names=FALSE,quote=FALSE)


##load the data of overlapped genes 21 genes
overlap <-read.delim("Tab3_overlap.txt", header = TRUE) # counts can be read into a data.frame

# view of data loaded
head(overlap)
dim(overlap)
colnames(overlap)


head(Tab1_Top.deg)
dim(Tab1_Top.deg)
colnames(overlap)


overlap <- as.data.frame(overlap)
head(overlap, n=10)
dim(overlap)


#merge the two tables (DEG gene from Lab strain with annoted table)
overlap <- merge(overlap,Tab1, by="GeneID")
head(overlap, n=10)
dim(overlap)


head(annot, n=10)
dim(annot)


#merge the two tables for LogFC, FDR, padj,.. with annoted table
overlap <- merge(overlap,annot, by="GeneID")
head(overlap, n=10)
dim(overlap)

#save the table in csv format
write.csv(overlap,"Tab3_overlap.deg.csv",row.names=FALSE,quote=FALSE)

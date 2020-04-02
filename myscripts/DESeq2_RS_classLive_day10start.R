## Set your working directory
setwd("C:/Users/kmb057/Documents/GitHub/EcologicalGenomics")

## Import the libraries that we're likely to need in this session
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)  ### First: BiocManager::install("vsn") AND BiocManager::install("hexbin")

## Import the counts matrix
countsTable <- read.table("C:/Users/kmb057/Documents/GitHub/EcologicalGenomics/RS_counts_samples/RS_cds2kb_countsMatrix.txt", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)
countsTableRound <- round(countsTable) # Need to round because DESeq wants only integers
head(countsTableRound)

## Import the samples description table - links each sample to factors of the experimental design.
# Need the colClasses otherwise imports "day" as numeric which DESeq doesn't like, coula altneratively change to d0, d5, d10
conds <- read.delim("C:/Users/kmb057/Documents/GitHub/EcologicalGenomics/RS_counts_samples/RS_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1, colClasses=c('factor', 'factor', 'factor', 'factor'))
head(conds)
dim(conds)

############ Try with only Day 10 data

# grep("10", names(countsTableRound), value = TRUE)
# day10countstable <- subset(countsTableRound, grep("10", names(countsTableRound), value = TRUE)) #doesn't work has to be logical

day10countstable <- countsTableRound %>% select(contains("10"))
dim(day10countstable)

conds10<- subset(conds, day=="10")
dim(conds10)
head(conds10)

## Let's see how many reads we have from each sample:
colSums(day10countstable)
mean(colSums(day10countstable))
barplot(colSums(day10countstable), las=3, cex.names=0.5,names.arg = substring(colnames(day10countstable),1,13))
abline(h=mean(colSums(day10countstable)), col="blue", lwd =2)

# What's the average number of counts per gene
rowSums(day10countstable)
mean(rowSums(day10countstable))
median(rowSums(day10countstable))
# wow! This shows dispersion across genes - differences in magnitude of expression

# What's the average number of counts per gene per sample
apply(countsTableRound,2,mean)

## Create a DESeq object and define the experimental design here with the tilde

dds <- DESeqDataSetFromMatrix(countData = day10countstable, colData = conds10, 
                              design = ~ climate + treatment + climate:treatment)
dim(dds)
# [1] 66408    30

# Filter out genes with few reads 

dds <- dds[rowSums(counts(dds)) > 30]
dim(dds)
# 24300    30
## Run the DESeq model to test for differential gene expression: 
# 1) estimate size factors (per sample), 2) estimate dispersion (per gene), 
# 3) run negative binomial glm
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)
# Running the model: design = ~ climate + treatment + climate:treatment

# [1] "Intercept"            "climate_HD_vs_CW"     "treatment_D_vs_C"    
# [4] "treatment_H_vs_C"     "climateHD.treatmentD" "climateHD.treatmentH"


# Order and list and summarize results from specific contrasts
# Here you set your adjusted p-value cutoff, can make summary tables of the number of genes differentially expressed (up- or down-regulated) for each contrast
res <- results(dds, name="treatment_D_vs_C", alpha = 0.05)
res <- res[order(res$padj),]
head(res)

# log2 fold change (MLE): treatment H vs C 
# Wald test p-value: treatment H vs C 
# DataFrame with 6 rows and 6 columns
# baseMean    log2FoldChange             lfcSE              stat
# <numeric>         <numeric>         <numeric>         <numeric>
#   MA_172878g0010   15.8548874481417  2.26899213338594  0.44065452286252  5.14914068882447
# MA_107783g0020    6.6082118492291 -1.96824414729957 0.394042285152064 -4.99500744327481
# MA_28973g0010    18.8813749792546  -1.9664671947209 0.412333404130031 -4.76911929769524
# MA_10434037g0010   5.611769238156   2.1853605976071  0.49671691935347  4.39960974240937
# MA_10426002g0010 10.8980752578363 -1.20767724886483 0.283132420279441 -4.26541491671245
# MA_10429525g0010 60.5937645508928  1.17170086164903 0.281505776107537  4.16226223792093
# pvalue                padj
# <numeric>           <numeric>
#   MA_172878g0010   2.61682549726074e-07 0.00249278796869058
# MA_107783g0020   5.88334982423675e-07 0.00280223952128396
# MA_28973g0010    1.85033065465192e-06 0.00587541660540472
# MA_10434037g0010  1.0844572516541e-05  0.0258263494481425
# MA_10426002g0010 1.99531043147766e-05  0.0357522692442608
# MA_10429525g0010 3.15110173906545e-05  0.0357522692442608

summary(res)
# out of 23887 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 16, 0.067%
# LFC < 0 (down)     : 3, 0.013%
# outliers [1]       : 61, 0.26%
# low counts [2]     : 14300, 60%

res_interClimTreat <- results(dds, name="climateHD.treatmentD", alpha=0.05)
res_interClimTreat <- res_interClimTreat[order(res_interClimTreat$padj),]
head(res_interClimTreat)
# log2 fold change (MLE): treatment D vs C 
# Wald test p-value: treatment D vs C 
# DataFrame with 6 rows and 6 columns
# baseMean   log2FoldChange             lfcSE             stat
# <numeric>        <numeric>         <numeric>        <numeric>
#   MA_10257300g0010 20.9979917001674  6.3160877571623 0.761778438070059 8.29123986911983
# MA_444738g0020   23.5872071084088 2.60491779951534 0.331247263623815 7.86396775332654
# MA_57964g0010    7.89927388331396 5.39652442842906 0.688793826922627 7.83474563432065
# MA_75192g0010    37.9183573851468 5.81210235837303 0.768762098604424 7.56033936756775
# MA_10428616g0010  35.758883777048 3.82582283371241 0.510641392538861 7.49219097709799
# MA_7017g0010     64.7924705055064 2.64439151272771 0.360360746749988 7.33817857959526
# pvalue                 padj
# <numeric>            <numeric>
#   MA_10257300g0010 1.12074011571854e-16 1.84462615646114e-12
# MA_444738g0020   3.72153387856494e-15 2.57743989393094e-11
# MA_57964g0010    4.69792799185419e-15 2.57743989393094e-11
# MA_75192g0010    4.02019076294682e-14 1.65420799418354e-10
# MA_10428616g0010 6.77332669682435e-14 2.22964368206064e-10
# MA_7017g0010     2.16520151829775e-13 5.93950863161046e-10

summary(res_interClimTreat)
# out of 23887 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 678, 2.8%
# LFC < 0 (down)     : 424, 1.8%
# outliers [1]       : 61, 0.26%
# low counts [2]     : 7367, 31%




# MA plot
plotMA(res_interClimTreat,ylim=c(-10,10))

# PCA
vsd <- vst(dds, blind=FALSE)

data <- plotPCA(vsd,intgroup=c("climate","treatment","pop"),returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

data$treatment <- factor(data$treatment, levels=c("C","H","D"), labels = c("C","H","D"))
data$day <- factor(data$day, levels=c("0","5","10"), labels = c("0","5","10"))


ggplot(data, aes(PC1, PC2, color=pop, shape=climate)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

# Counts of specific top gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds, gene="MA_7017g0010", intgroup = (c("treatment","climate","day")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=climate)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size=3) +
  scale_x_discrete(limits=c("C","H","D"))
p

p <-ggplot(d, aes(x=treatment, y=count, color=climate)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p

# Heatmap of top 20 genes sorted by pvalue

library(pheatmap)
topgenes <- head(rownames(res_treatCD),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("treatment","climate","day")])
pheatmap(mat, annotation_col=df)






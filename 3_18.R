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

## Let's see how many reads we have from each sample:
colSums(countsTableRound)
mean(colSums(countsTableRound))
barplot(colSums(countsTableRound), las=3, cex.names=0.5,names.arg = substring(colnames(countsTableRound),1,13))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd =2)

#What is the average number of counts per gene?
rowSums(countsTableRound)
mean(rowSums(countsTableRound))
median(rowSums(countsTableRound))
#wow this shows dispersion across genes. Mean thrown off by few highly expressed genes. Differnces in magnitude of expression.


#What's the average number of counts per gene per sample?
apply(countsTableRound,2,mean)

## Create a DESeq object and define the experimental design here with the tilde
dds <- DESeqDataSetFromMatrix(countData = countsTableRound,colData = conds, design = ~ pop + day + treatment)
dim(dds)
#[1] 66408    76

# Filter out genes with few reads
dds <- dds[rowSums(counts(dds)) > 760]
dim(dds)
#[1] 23887    76 This is filtering to sum of 76 reads across all samples
#[1] 7884   76 This is filtering to sum of 760 reads across all samples


## Run the DESeq model to test for differential gene expression: 1) estimate size factors (per sample), 2) estimate dispersion (per gene), 3) run negative binomial glm

dds <- DESeq(dds)
resultsNames(dds)
#[1] "Intercept"           
#[2] "pop_BRU_05_vs_ASC_06"
#[3] "pop_CAM_02_vs_ASC_06"
#[4] "pop_ESC_01_vs_ASC_06"
#[5] "pop_JAY_02_vs_ASC_06"
#[6] "pop_KAN_04_vs_ASC_06"
#[7] "pop_LOL_02_vs_ASC_06"
#[8] "pop_MMF_13_vs_ASC_06"
#[9] "pop_NOR_02_vs_ASC_06"
#[10] "pop_XBM_07_vs_ASC_06"
#[11] "day_10_vs_0"         
#[12] "day_5_vs_0"          
#[13] "treatment_D_vs_C"    
#[14] "treatment_H_vs_C"


dds <- DESeqDataSetFromMatrix(countData = countsTableRound,colData = conds, design = ~ climate + day + treatment)
dim(dds)
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)
#[1] "Intercept"       
#[2] "climate_HD_vs_CW"
#[3] "day_10_vs_0"     
#[4] "day_5_vs_0"      
#[5] "treatment_D_vs_C"
#[6] "treatment_H_vs_C"


# Order and list and summarize results from specific contrasts
# Here you set your adjusted p-value cutoff, can make summary tables of the number of genes differentially expressed (up- or down-regulated) for each contrast
res <- results(dds,alpha = 0.05)
res <- res[order(res$padj)]
head(res)
#log2 fold change (MLE): treatment H vs C 
#Wald test p-value: treatment H vs C 
#DataFrame with 6 rows and 6 columns
#baseMean
#<numeric>
#   MA_10000213g0010                   0
# MA_10000405g0010 0.00927779777599038
# MA_10000516g0010   0.306727408755282
# MA_10001015g0010  0.0953769494456185
# MA_10001337g0010    10.5864076666598
# MA_10002583g0010  0.0720400016104675
# log2FoldChange
# <numeric>
#   MA_10000213g0010                  NA
# MA_10000405g0010  -0.678234333988712
# MA_10000516g0010   0.674578980909189
# MA_10001015g0010 -0.0643026331356064
# MA_10001337g0010  -0.176203680815192
# MA_10002583g0010  -0.112494012427474
# lfcSE
# <numeric>
#   MA_10000213g0010                NA
# MA_10000405g0010  3.56904296103673
# MA_10000516g0010 0.989443886861717
# MA_10001015g0010  3.56705966952029
# MA_10001337g0010 0.426103598391028
# MA_10002583g0010  3.56871636664819
# stat
# <numeric>
#   MA_10000213g0010                  NA
# MA_10000405g0010  -0.190032549731959
# MA_10000516g0010   0.681775884278587
# MA_10001015g0010 -0.0180267893147562
# MA_10001337g0010  -0.413523099735695
# MA_10002583g0010 -0.0315222620320287
# pvalue
# <numeric>
#   MA_10000213g0010                NA
# MA_10000405g0010 0.849283624250601
# MA_10000516g0010 0.495380675441948
# MA_10001015g0010 0.985617482098602
# MA_10001337g0010 0.679223401788245
# MA_10002583g0010 0.974853038430995
# padj
# <numeric>
#   MA_10000213g0010                NA
# MA_10000405g0010                NA
# MA_10000516g0010                NA
# MA_10001015g0010                NA
# MA_10001337g0010 0.933204889831224
# MA_10002583g0010

summary(res)
# out of 53066 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 15, 0.028%
# LFC < 0 (down)     : 6, 0.011%
# outliers [1]       : 307, 0.58%
# low counts [2]     : 45381, 86%
# (mean count < 8)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

res_treatCD <- results(dds, name="treatment_D_vs_C",alpha=0.05)
res_treatCD <- res_treatCD[order(res_treatCD$padj),]
head(res_treatCD)
# log2 fold change (MLE): treatment D vs C 
# Wald test p-value: treatment D vs C 
# DataFrame with 6 rows and 6 columns
# baseMean
# <numeric>
#   MA_10257300g0010 20.9979917001674
# MA_444738g0020   23.5872071084088
# MA_57964g0010    7.89927388331396
# MA_75192g0010    37.9183573851468
# MA_10428616g0010  35.758883777048
# MA_7017g0010     64.7924705055064
# log2FoldChange
# <numeric>
#   MA_10257300g0010  6.3160877571623
# MA_444738g0020   2.60475479951658
# MA_57964g0010    5.39739873586083
# MA_75192g0010    5.81210235837303
# MA_10428616g0010 3.82582767481363
# MA_7017g0010     2.64414218710183
# lfcSE
# <numeric>
#   MA_10257300g0010 0.761778438070059
# MA_444738g0020   0.330944830911607
# MA_57964g0010    0.694994697724841
# MA_75192g0010    0.768762098604424
# MA_10428616g0010 0.510533466865097
# MA_7017g0010     0.358927431461684
# stat
# <numeric>
#   MA_10257300g0010 8.29123986911983
# MA_444738g0020   7.87066168201399
# MA_57964g0010    7.76610059548648
# MA_75192g0010    7.56033936756775
# MA_10428616g0010 7.49378429254779
# MA_7017g0010     7.36678769949098
# pvalue
# <numeric>
#   MA_10257300g0010 1.12074011571857e-16
# MA_444738g0020   3.52770431249132e-15
# MA_57964g0010    8.09392779867025e-15
# MA_75192g0010    4.02019076294688e-14
# MA_10428616g0010 6.69156646447227e-14
# MA_7017g0010     1.74788498523708e-13
# padj
# <numeric>
#   MA_10257300g0010 1.84910911692407e-12
# MA_444738g0020   2.91017967258971e-11
# MA_57964g0010    4.45139049167535e-11
# MA_75192g0010    1.65822818494651e-10
# MA_10428616g0010 2.20808310194656e-10
# MA_7017g0010     4.80639239523776e-10

summary(res_treatCD)
# out of 53066 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 667, 1.3%
# LFC < 0 (down)     : 414, 0.78%
# outliers [1]       : 307, 0.58%
# low counts [2]     : 36260, 68%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

##### Data visualization #####
# MA plot
plotMA(res_treatCD,ylim=c(-3,3))

# PCA
vsd <- vst(dds,blind=FALSE)

data <- plotPCA(vsd,intgroup=c("climate","treatment","day"),returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))

ggplot(data, aes(PC1, PC2, color=day, shape=treatment)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

ggplot(data, aes(PC1, PC2, color=climate, shape=treatment)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

# Counts of specific top gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds, gene="MA_10426407g0030", intgroup = (c("treatment","climate")), returnData=TRUE)
d



p <-ggplot(d, aes(x=climate, y=count, shape=climate, colour = treatment)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size=3) +
  scale_x_discrete(limits=c("CW","HD"))
p



p <-ggplot(d, aes(x=treatment, y=count, shape=climate)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p

# Heatmap of top 20 genes sorted by pvalue


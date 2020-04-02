library(methylKit)
library(tidyverse)
library(ggplot2)
library(pheatmap)

# first, we want to read in the raw methylation calls with methylkit

# set directory with absolute path (why is this necessary? I have no idea, but gz files wont work with relative paths)
dir <- "C:/Users/kmb057/Documents/EcologicalGenomics/Epigenetics_data/Epigenetics_data"
# read in the sample ids
samples <- read.table("C:/Users/kmb057/Documents/EcologicalGenomics/Epigenetics_data/Epigenetics_data/sample_id.txt",header=FALSE)
# now point to coverage files
files <- file.path(dir,samples$V1)
all(file.exists(files))
# convert to list
file.list <- as.list(files)
# get the names only for naming our samples
nmlist <- as.list(gsub("_1_bismarck_bt2_pe.bismarck.cov.gz","",samples$V1))
# use methRead to read in the coverage files
myobj <- methRead(location= file.list,
                  sample.id =   nmlist,
                  assembly = "atonsa", # this is just a string. no actual database
                  dbtype = "tabix",
                  context = "CpG",
                  resolution = "base",
                  mincov = 20,
                  treatment = 
                    c(0,0,0,0,
                      1,1,1,1,
                      2,2,2,2,
                      3,3,3,3,
                      4,4,4,4),
                  pipeline = "bismarkCoverage")
######
# visualize coverage and filter
######

# We can look at the coverage for individual samples with getCoverageStats()
getCoverageStats(myobj[[1]],plot=TRUE)

# and can plot all of our samples at once to compare.

# filter samples by depth with filterByCoverage()
filtered.myobj <- filterByCoverage(myobj,
                                lo.count=20,lo.perc=NULL,
                                hi.count=NULL,hi.perc=97.5,
                                db.dir="C:/Users/kmb057/Documents/EcologicalGenomics/Epigenetics_data/Epigenetics_data")

######
# merge samples
######

#Note! This takes a while and we're skipping it

# use unite() to merge all the samples. We will require sites to be present in each sample or else will drop it

meth <- unite(filtered.myobj,mc.cores=3,suffix="united",
              db.dir="C:/Users/kmb057/Documents/EcologicalGenomics/Epigenetics_data/Epigenetics_data")

meth <- methylKit:::readMethylBaseDB(
  dbpath = "C:/Users/kmb057/Documents/EcologicalGenomics/Epigenetics_data/Epigenetics_data/methylBase_united.txt.bgz",
  dbtype = "tabix",
  sample.id =   unlist(nmlist),
  assembly = "atonsa", # this is just a string. no actual database
  context = "CpG",
  resolution = "base",
  treatment = c(0,0,0,0,
                1,1,1,1,
                2,2,2,2,
                3,3,3,3,
                4,4,4,4),
  destrand = FALSE)

# percMethylation() calculates the percent methylation for each site and sample
pm <- percMethylation(meth)
#plot methylation histograms
ggplot(gather(as.data.frame(pm)),aes(value))+
  geom_histogram(bins=10,color="black",fill="grey")+
  facet_wrap(~key)

# calculate and plot mean methylation
sp.means <- colMeans(pm)

p.df <- data.frame(sample=names(sp.means),
                   group=substr(names(sp.means), 1,6),
                   methylation= sp.means)

ggplot(p.df,aes(x=group,y=methylation,color=group))+
  stat_summary(color="black")+geom_jitter(width=0.1,size=3)

clusterSamples(meth, dist="correlation", method="ward.D", plot=TRUE)

PCASamples(meth, screeplot=TRUE)
PCASamples(meth, screeplot=FALSE)

# subset with reorganize()

meth_sub <- reorganize(meth,  sample.ids= (c("AA_F00_1_1_bismark_bt2_pe.bismark.cov.gz","AA_F00_2_1_bismark_bt2_pe.bismark.cov.gz","AA_F00_3_1_bismark_bt2_pe.bismark.cov.gz", "AA_F00_4_1_bismark_bt2_pe.bismark.cov.gz",
                                             "AA_F25_3_1_bismark_bt2_pe.bismark.cov.gz","AA_F25_2_1_bismark_bt2_pe.bismark.cov.gz","AA_F25_3_1_bismark_bt2_pe.bismark.cov.gz","AA_F25_4_1_bismark_bt2_pe.bismark.cov.gz")), 
                       treatment=c(0,0,0,0,1,1,1,1),
                       save.db=FALSE)

# calculate differential methylation

myDiff=calculateDiffMeth(meth_sub,
                         overdispersion="MN",
                         mc.cores=1,
                         suffix = "full_model", adjust="qvalue",test="Chisq")

# where MN corrects for overdispersion
# fit a logistic regression to methylation values where explanatory variable is the treatment (case or control). 
# and we compare the fit of the model with explanatory variable vs the null model (no explanatory variable) 
#and ask if the fit is better using a Chisq test. 
# the methylation proportions are weighted by their coverage, as in a typical logistic regression. Note that in theory you could enter these as two column success and failure data frame, which is common in logistic regressions.

# use overdispersion: Chisq without overdispersion finds more true positives, but many more false positives. good compromise is overdispersion with Chisq. reduced true pos, but really reduces false pos rate.

# get all differentially methylated bases
myDiff=getMethylDiff(myDiff,difference=10,qvalue=0.05)

# we can visualize the changes in methylation frequencies quickly.
hist(getData(myDiff)$meth.diff)

# get hyper methylated bases
hyper=getMethylDiff(myDiff,difference=10,qvalue=0.05,type="hyper")
#
# get hypo methylated bases
hypo=getMethylDiff(myDiff,difference=10,qvalue=0.05,type="hypo")

#heatmaps first

# get percent methylation matrix
pm <- percMethylation(meth_sub)

# make a dataframe with snp id's, methylation, etc.
sig.in <- as.numeric(row.names(myDiff))
pm.sig <- pm[sig.in,]

# add snp, chr, start, stop

din <- getData(myDiff)[,1:3]
df.out <- cbind(paste(getData(myDiff)$chr, getData(myDiff)$start, sep=":"), din, pm.sig)
colnames(df.out) <- c("snp", colnames(din), colnames(df.out[5:ncol(df.out)]))
df.out <- (cbind(df.out,getData(myDiff)[,5:7]))

####
# heatmap
####

my_heatmap <- pheatmap(pm.sig,
                       show_rownames = FALSE)

ctrmean <- rowMeans(pm.sig[,1:4])

h.norm <- (pm.sig-ctrmean)

my_heatmap <- pheatmap(h.norm,
                       show_rownames = FALSE)

##### if you want to change colors. only because I don't love the default colors.

paletteLength <- 50
myColor <- colorRampPalette(c("cyan1", "black", "yellow1"))(paletteLength)
myBreaks <- c(seq(min(h.norm), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(h.norm)/paletteLength, max(h.norm), length.out=floor(paletteLength/2)))

my_heatmap <- pheatmap(h.norm,
                       color=myColor, 
                       breaks=myBreaks,
                       show_rownames = FALSE)

#####
#let's look at methylation of specific gene or snp
####

df.out
df.plot <- df.out[,c(1,5:12)] %>% pivot_longer(-snp, values_to = "methylation")
df.plot$group <- substr(df.plot$name,1,2)
head(df.plot)

# looking at snp LS049205.1:248
# if you choose a different snp, you can create different plots.

df.plot %>% filter(snp=="LS041600.1:288") %>% 
  ggplot(., aes(x=group, y=methylation, color=group, fill=group)) +
  stat_summary(fun.data = "mean_se", size = 2) +
  geom_jitter(width = 0.1, size=3, pch=21, color="black")

## write bed file for intersection with genome annotation

write.table(file = "C:/Users/kmb057/Documents/EcologicalGenomics/Epigenetics_data/Epigenetics_data/diffmeth.bed",
            data.frame(chr= df.out$chr, start = df.out$start, end = df.out$end),
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")




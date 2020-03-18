setwd("~/Documents/GitHub/EcologicalGenomics/myresults")

list.files()

SFS <- scan("XFS_outFold.sfs")

sumSFS <- sum(SFS)


pctPoly = 100*(1-(SFS[1]/sumSFS))

plotSFS <- SFS[-c(1,length(SFS))]

barplot(plotSFS)

div <- read.table("XFS_folded_allsites.thetas.idx.pestPG")

colnames(div)=c("window","chrname","wincenter","tW","tP","tF","tH","tL","tajD","fulif","fuliD","fayH","zengsE","numSites")

div$tWpersite=div$tW/div$numSites
div$tPpersite=div$tP/div$numSites

hist(div$tPpersite-div$tWpersite)

pdf("XFS_diversity_stats.pdf")
par(mfrow=c(2,2))
hist(div$tWpersite, col="gray",xlab="Theta-W",main="")
hist(div$tPpersite, col="gray",xlab="Theta-Pi",main="")
hist(div$tajD, col="gray",xlab="Tajima's D",main="")
barplot(plotSFS)
dev.off()

summary(div)

head(div)

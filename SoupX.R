###### Perform SoupX  
# ZY
# 3 Nov 2023
# refer to: https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html

library(SoupX)
library(Seurat)
library(Drop)

### quick methods-----
sc1 <- load10X("/data/share_for_user/1207_zhangjunjLab/OHSV1/outs")
sc2 <- load10X("/data/share_for_user/1207_zhangjunjLab/OHSV2/outs")
sc3 <- load10X("/data/share_for_user/1207_zhangjunjLab/OHSV3/outs")
sc4 <- load10X("/data/share_for_user/1207_zhangjunjLab/PBS1/outs")
sc5 <- load10X("/data/share_for_user/1207_zhangjunjLab/PBS2/outs")
sc6 <- load10X("/data/share_for_user/1207_zhangjunjLab/PBS3/outs")

png("/home/zhuyong/Wangshaowei_10X_6samples/results/OHSV1_rho.png")
sc1 <- autoEstCont(sc1)
dev.off()

png("/home/zhuyong/Wangshaowei_10X_6samples/results/OHSV2_rho.png")
sc2 <- autoEstCont(sc2)
dev.off()

png("/home/zhuyong/Wangshaowei_10X_6samples/results/OHSV3_rho.png")
sc3 <- autoEstCont(sc3)
dev.off()

png("/home/zhuyong/Wangshaowei_10X_6samples/results/PBS1_rho.png")
sc4 <- autoEstCont(sc4)
dev.off()

png("/home/zhuyong/Wangshaowei_10X_6samples/results/PBS2_rho.png")
sc5 <- autoEstCont(sc5)
dev.off()

png("/home/zhuyong/Wangshaowei_10X_6samples/results/PBS3_rho.png")
sc6 <- autoEstCont(sc6)
dev.off()

out1 <-  adjustCounts(sc1)
DropletUtils:::write10xCounts("/data/share_for_user/1207_zhangjunjLab/OHSV1/outs/counts_soupX", out1)
out2 <-  adjustCounts(sc2)
DropletUtils:::write10xCounts("/data/share_for_user/1207_zhangjunjLab/OHSV2/outs/counts_soupX", out2)
out3 <-  adjustCounts(sc3)
DropletUtils:::write10xCounts("/data/share_for_user/1207_zhangjunjLab/OHSV3/outs/counts_soupX", out3)
out4 <-  adjustCounts(sc4)
DropletUtils:::write10xCounts("/data/share_for_user/1207_zhangjunjLab/PBS1/outs/counts_soupX", out4)
out5 <-  adjustCounts(sc5)
DropletUtils:::write10xCounts("/data/share_for_user/1207_zhangjunjLab/PBS2/outs/counts_soupX", out5)
out6 <-  adjustCounts(sc6)
DropletUtils:::write10xCounts("/data/share_for_user/1207_zhangjunjLab/PBS3/outs/counts_soupX", out6)

### step by step methods(-ing)----
# add extra meta data to the SoupChannel object
sc <- setClusters(sc, setNames(metaData$Cluster, rownames(metaData)))
sc <- setDR(sc, metaData[colnames(sc$toc), c("RD1", "RD2")])

# Visual sanity checks
dd <-  metaData[colnames(sc$toc), ]
mids <-  aggregate(cbind(RD1, RD2) ~ Annotation, data = dd, FUN = mean)
gg <-  ggplot(dd, aes(RD1, RD2)) + geom_point(aes(colour = Annotation), size = 0.2) + 
  geom_label(data = mids, aes(label = Annotation)) + ggtitle("Annotation") + 
  guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

dd$geneA <-  sc$toc["geneA", ]
gg <-  ggplot(dd, aes(RD1, RD2)) + geom_point(aes(colour = geneA > 0))
plot(gg)

gg <-  plotMarkerMap(sc, "geneA")
plot(gg)

# estimate the contamination fraction
sc <- setContaminationFraction(sc, 0.2) ###set the contamination fraction to 20% for all cells.

#pick soup specific genes manually
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
plotMarkerDistribution(sc)





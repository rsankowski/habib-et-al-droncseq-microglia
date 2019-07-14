#RaceID4 
library(tidyverse)
library(viridis)
library(RaceID)
library(Matrix)

date = Sys.Date()
#download counts file from url: https://stackoverflow.com/questions/28986150/downloading-and-extracting-gz-data-file-using-r
#dronc dataset

url <- "https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz"
tmp <- tempfile()
##
download.file(url,tmp)
dronc <- read.delim(gzfile(tmp), sep = '\t', stringsAsFactors = F, skip = 2)
genes <- dronc$Description
dronc <- Matrix::Matrix(as.matrix(dronc[,3:ncol(dronc)]), sparse = T)
rownames(dronc) <- genes
which(genes == "TMEM119")
hist(dronc[which(genes == "TMEM119"),])
hist(dronc[which(genes == "SLC2A5"),])
hist(dronc[which(genes == "P2RY12"),])
sum(dronc[which(genes == "TMEM119"),]>5000)
sum(dronc[which(genes == "SLC2A5"),]>5000)
sum(dronc[which(genes == "SLC2A5"),]>2000)

micr <- dronc[,which(dronc[which(genes == "TMEM119"),]>750 & dronc[which(genes == "SLC2A5"),]>750 & dronc[which(genes == "SLC2A5"),]>750)]

save(micr, file = "data/Habib-et-al-droncseq-microglia.Robj")

load("data/Habib-et-al-droncseq-microglia.Robj")

prdata <- micr

sc <- SCseq(prdata)
# filtering of expression data
a <- apply(prdata, 2, sum)
sc <- filterdata(sc, mintotal=quantile(a, 0.1)) # exlcude the lower quartile of the cells
sc <- CCcorrect(sc, 
                dimR = T, 
                nComp = 20,
                CGenes = c('JUN',
                           'FOS',
                           'ZFP36',
                           'HSPA1A|HSPA1B',
                           'DUSP1',
                           'EGR1',
                           'MALAT1'))

sc <- compdist(sc,metric="pearson")
sc <- clustexp(sc) 

plotsaturation(sc,disp=FALSE)
plotsaturation(sc,disp=TRUE)
plotjaccard(sc)

sc <- clustexp(sc,cln=16,sat=FALSE) 
sc <- findoutliers(sc)
plotbackground(sc)
plotsensitivity(sc)
plotoutlierprobs(sc)
ord_clust <- clustheatmap(sc)
save(ord_clust, file = 'data/ord_clust.Robj')

pdf(paste0('plots/heatmaps/clustheatmap.pdf'))
clustheatmap(sc, final = T)
dev.off()

sc <- comptsne(sc)
sc <- compfr(sc,knn=10)

plotmap(sc)
plotmap(sc,fr=TRUE)
dev.off()

name2id <- function(x,id) {
  ##  id[sub("\\_\\_chr\\w+","",id) %in% x]
  n <- c()
  for ( j in x ){ n <- append(n,id[grep(paste(j,"(\\_\\_chr\\w+|$|\\|)",sep=""),id)])
  }
  n
}

plotexpmap(sc,name2id("MRC1", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("LYVE1", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("CD163", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Tmem119", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("CX3CR1", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("PTPRC", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("CD3E", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("ITGAM", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("CD8A", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("CD4", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("P2RY12", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("SLC2A5", rownames(sc@ndata)),logsc=F,fr=F)

dg <- clustdiffgenes(sc,4,pvalue=.01)
head(dg,25)
types <- sub("(\\_\\d+)$","", colnames(sc@ndata))
genes <- head(rownames(dg)[dg$fc>1],10)
plotmarkergenes(sc,genes,samples=types)

#Save sc file
save(sc, file = 'data/sc.Robj')

micr_ids <- names(sc@cpart)[sc@cpart %in% c(17,9,19,10)]
write_csv(as.data.frame(micr_ids), "microglia-cell-ids.csv")

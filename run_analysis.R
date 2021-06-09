### load libraries
library(vroom)
library(tibble)
library(Rsamtools)
library(GenomicRanges)
library(DESeq2)
library(lme4)
library(lmerTest)
library(multcomp)
library(ggplot2)
library(mda)
library(gplots)

## load regression custom functions
source("src/functions.R")


## define directory paths
mainDir <- $PWD
dataDir <- paste(mainDir, 'data/raw/', sep='')
dataProcDir <- paste(mainDir, 'data/processed/', sep='')
repDir <- paste(mainDir, 'output/', sep='')
tabDir <- paste(repDir, 'tables/', sep='')
figDir <- paste(repDir, 'figs/', sep='')
# temporary directory
tmpDir <- "/local_scratch/malvarez/chrysper/analysis/"



### load info on samples
sampleinfo <- read.table(paste0(dataDir, "sampleinfo.txt"), sep="\t", header=T, stringsAsFactors=F)[,1:10]
sampleinfo$filename <- paste0(dataDir, "fastas/180329/0", sampleinfo$FGC.Pool.ID, "_", sampleinfo$CRG.Sample.ID, "_", substr(sampleinfo$barcode.seq, 1, 7), ".fastq.gz")
colnames(sampleinfo)[colnames(sampleinfo) == "group"] <- "sampleName"
sampleinfo$group <- paste0(c(rep("tumor", 12), "tcbt", "tcbt", "cu", "cdt", "cdr.res"), "_", gsub(" ", "", sampleinfo$treatment))


### load info on guides
f1 <- paste0(dataDir, "broadgpp-brunello-library-contents.txt")
gl <- read.table(f1, stringsAsFactors=F, sep="\t", header=T)
rownames(gl) <- paste0(gl[,2], "-", 1:nrow(gl))


### load matrix of raw read counts
mat = vroom(paste0(dataDir, "mat.counts.csv")) %>%
  column_to_rownames("sgRNA") %>%
  # all counts are 0 on column cdt, dropped from the analysis
  select(-cdt) %>%
  as.matrix()


### normalize count matrix (3 ways)

## rlog
rlog = rlog(mat)

## varianceStabilizingTransformation
vst = varianceStabilizingTransformation(mat)

## log2-normalized-reads-per-million for each sgRNA
norm = log2( ( t( t(mat) / colSums(mat) ) * 1e6 ) + 1)
# another version (Camille's)
norm2 = log2( t( t(mat + 1) / colSums(mat + 1) ) * 1e7 )

# print normalized sample distributions
pdf(paste0(figDir, "normalization.pdf"), width=15, height=6)
par(mfrow=c(1,3))
boxplot(vst, outline=F, las=2, main="vst")
boxplot(rlog, outline=F, las=2, main="rlog")
boxplot(norm, outline=F, las=2, main="shift log")
boxplot(norm2, outline=F, las=2, main="shift log Camille")
dev.off()


## create Biobase's 'ExpressionSets' (esets)
esets <- list()

eset <- new("ExpressionSet", exprs=rlog)
pData(eset) <- sampleinfo[match(colnames(eset), sampleinfo$sampleName),]
eset$group <- factor(eset$group)
rownames(eset) <- rownames(mat)
fData(eset) <- gl[rownames(eset), ]
esets[["rlog"]] <- eset

eset <- new("ExpressionSet", exprs=vst)
pData(eset) <- sampleinfo[match(colnames(eset), sampleinfo$sampleName),]
eset$group <- factor(eset$group)
rownames(eset) <- rownames(mat)
fData(eset) <- gl[rownames(eset), ]
esets[["vst"]] <- eset

eset <- new("ExpressionSet", exprs=norm)
pData(eset) <- sampleinfo[match(colnames(eset), sampleinfo$sampleName),]
eset$group <- factor(eset$group)
rownames(eset) <- rownames(mat)
fData(eset) <- gl[rownames(eset), ]
esets[["norm"]] <- eset

eset <- new("ExpressionSet", exprs=norm2)
pData(eset) <- sampleinfo[match(colnames(eset), sampleinfo$sampleName),]
eset$group <- factor(eset$group)
rownames(eset) <- rownames(mat)
fData(eset) <- gl[rownames(eset), ]
esets[["norm2"]] <- eset

save(esets, file=paste0(dataDir, "esets.RData"))
#load(file=paste0(dataDir, "esets.RData"))


## PCAs
for(i in names(esets)){
  pr <- prcomp(t(exprs(esets[[i]])))
  pcvar <- round(pr$sdev^2/sum(pr$sdev^2)*100, 2)
  ##
  pdf(paste0(figDir, "pca.", i, ".pdf"), width=9, height=9)
  par(mfrow=c(3, 3))
  for(i in 2:3){
    plot(pr$x[,1], pr$x[,i], col=rainbow(5)[esets[[i]]$group], xlab=paste0("PC1 (", pcvar[1], "%)"), ylab=paste0("PC", i, " (", pcvar[i], "%)"))
    text(pr$x[,1], pr$x[,i], labels=esets[[i]]$sampleName, pos=4, cex=0.6, col=rainbow(5)[esets[[i]]$group], srt=90)
  }
  plot(pr$x[,1], apply(exprs(esets[[i]]), 2, median), col=rainbow(5)[esets[[i]]$group], xlab=paste0("PC1 (", pcvar[1], "%)"), ylab="Col medians")
  plot(pr$x[,1], colSums(mat), col=rainbow(5)[esets[[i]]$group], xlab=paste0("PC1 (", pcvar[1], "%)"), ylab="Col sums")
  plot(pr$x[,2], colSums(mat), col=rainbow(5)[esets[[i]]$group], xlab=paste0("PC2 (", pcvar[2], "%)"), ylab="Col sums")
  plot(pr$x[,3], colSums(mat), col=rainbow(5)[esets[[i]]$group], xlab=paste0("PC3 (", pcvar[3], "%)"), ylab="Col sums")
  plot(pr$x[,1], colMeans(mat == 0), col=rainbow(5)[esets[[i]]$group], xlab=paste0("PC1 (", pcvar[1], "%)"), ylab="% zeros")
  plot(pr$x[,2], colMeans(mat), col=rainbow(5)[esets[[i]]$group], xlab=paste0("PC2 (", pcvar[2], "%)"), ylab="% zeros")
  plot(pr$x[,3], colMeans(mat), col=rainbow(5)[esets[[i]]$group], xlab=paste0("PC3 (", pcvar[3], "%)"), ylab="% zeros")
  dev.off()
}



### differential lethality

### Work with VST-normalized counts
eset <- esets[["vst"]]
sampleNames(eset) <- eset$sampleName
eset$mouse[13:16] <- 10:13
colnames(fData(eset))[2] <- "gene"
eset <- eset[!is.na(fData(eset)$gene),]
eset$treatment[13:14] <- "NA"; eset$treatment[16] <- "drug treatment"; eset$treatment <- as.factor(eset$treatment)

conts =  mcp(group = c("tumor_untreated - tcbt_NA = 0",
                       "tumor_drugtreatment - tcbt_NA = 0",
                       "tumor_drugtreatment - tumor_untreated = 0"))


## run regression (loop across gene names)
allres <- mclapply(unique(fData(eset)$gene), function(g){
  df <- buildDF(g, eset, ref="tcbt_NA", cond="group")
  y <- fitModel(df, conts, cond="group")
  m1 <- tapply(df$counts, df$group, mean)
  list(df=df, res=y, mns=m1)
}, mc.cores=30)
names(allres) <- unique(fData(eset)$gene)

save(allres, file=paste0(dataDir, "allres.RData"))


## parse results
cres <- lapply(1:length(conts$group), function(i){
  y <- do.call(rbind, lapply(allres, function(x) as.matrix(unname(x$res[i,, drop=F]))))
  colnames(y) <- colnames(allres[[1]]$res)
  mns <- do.call(rbind, lapply(allres, "[[", "mns"))
  y <- cbind(y, mns)
  rownames(y) <- names(allres)
  y
})
names(cres) <- rownames(allres[[1]]$res)

cres <- lapply(cres, function(x) {
  x <- data.frame(x)
  x$pval.adj <- p.adjust(x$pvals, method="BH")
  colnames(x) <- sub("ci.", "", colnames(x))
  x
})

mns <- cres[[1]][,5:9]

allres <- cbind(mns, untrVsTCBT = cres[[1]][, c(1:4, 10)], treatVsTCBT = cres[[2]][, c(1:4, 10)], treatVsuntr = cres[[3]][, c(1:4, 10)])

# write tables
write.table(allres, file=paste0(tabDir, "results.xls"), quote=F, row.names=T, sep="\t")

for(i in names(cres)) write.table(cres[[i]], file=paste0(tabDir, i, ".results.xls"), quote=F, row.names=T, sep="\t")


## detect top 50 genes -sorted by max(|estimate|)- with significant regression (|estimate|>4 and FDR<0.05)
sel <- unique(unlist(lapply(cres, function(x){
  x <- x[order(abs(x$Estimate), decreasing=T),]
  x <- x[x$pval.adj < 0.05 & abs(x$Estimate) > 4 ,]
  rownames(x)[1:min(nrow(x), 50)]
})))

mycols <- c(1, '#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506')
mycols5 <- mycols[c(1, 3, 5, 7, 9)]

system(paste0("mkdir ", figDir, "Stripcharts"))


## plots for the selected genes

# 2 barplots per gene
for(i in sel){
  
  y <- buildDF(i, eset, ref="tcbt_NA")
  
  pdf(paste0(figDir, "Stripcharts/", i, ".grna.pdf"), width=10, height=8)
  q = ggplot(data = y,
             aes(x = group,
                 y = counts,
                 colour = gRNA)) +
    geom_boxplot() +
    scale_colour_manual(values=mycols5) +
    ggtitle(i)
  print(q)
  dev.off()
  
  pdf(paste0(figDir, "Stripcharts/", i, ".group.pdf"), width=10, height=8)
  q = ggplot(data = y,
             aes(x = group,
                 y = counts,
                 colour = group)) +
    geom_boxplot() +
    scale_colour_manual(values=mycols5)
  print(q)
  dev.off()
}

# heatmap
heatmap.2(exprs(eset[fData(eset)$gene %in% sel,]), trace="none", scale="row", ColSideColors=c("red", "orange", "gray")[factor(eset$treat)])

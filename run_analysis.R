
### load libraries
library(magrittr)
library(vroom)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(writexl)
library(ggplot2)
library(MASS)
# library(Rsamtools)
# library(GenomicRanges)
library(DESeq2)
# library(multcomp)
# library(mda)
# library(gplots)
library(conflicted)

# avoid function conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")

## load regression custom functions
source("src/functions.R")



## define and create directories

# project directory
mainDir = getwd()

# processed data directory
dataDir = paste0(mainDir, 'data/processed/') ; dir.create(dataDir, recursive = TRUE)

# output directories
outDir = paste0(mainDir, 'output/')
tabDir = paste0(outDir, 'tables/') ; dir.create(tabDir, recursive = TRUE)
figDir = paste0(outDir, 'figs/') ; dir.create(figDir, recursive = TRUE)

# temporary directory
tmpDir = "/local_scratch/malvarez/chrysper/tmp/"



#### load data

# raw data's root path
rawPath = "/g/strcombio/fsupek_data/CRISPR/"


### load tables of raw read counts (3 different experiments)

## our APOBEC3A experiment
# cell lines --> A549, LXF289, H358
# treatments --> DOX (A3A active) vs. control -- some also have ATRi
# altered genotypes --> TP53-KO, TP53-mut, HMCES-KO
A3A_Dir = paste0(rawPath, "2_processed_data/apobec/MAGeCK_MLE_analyses/raw_counts_combined/")
A3A_raw_counts = vroom(paste0(A3A_Dir, "A549_LXF289_H358__raw_counts.tsv")) %>%
  as.matrix()
# load sample info 
A3A_sampleinfo = vroom(paste0(A3A_Dir, "sampleinfo.tsv"))
# more info
#/g/strcombio/fsupek_data/CRISPR/1_raw_data/0_INFO/paths_info_A549_LXF289_H358.xlsx

## Travis' ATRi/TMZ experiment
# cell lines --> SW480, SW620, HT29
# treatments --> ATRi, TMZ vs. control
# altered genotypes --> ARID1A-KO, MSH6-KO
ATRi_TMZ_Dir = paste0(rawPath, "/g/strcombio/fsupek_data/CRISPR/3_reports/stracker_lab/isabel_ATRi_TMZ/raw_counts_combined/")
ATRi_TMZ_raw_counts = vroom(paste0(ATRi_TMZ_Dir, "SW480_SW620_HT29__ARID1A_MSH6__raw_counts.tsv")) %>%
  as.matrix()
# load sample info 
ATRi_TMZ_sampleinfo = vroom(paste0(ATRi_TMZ_Dir, "sampleinfo.tsv"))
# more info
#readLines("/g/strcombio/fsupek_data/CRISPR/1_raw_data/0_INFO/paths_info_isabel.txt")

## Jain et al oxygen experiment
# cell lines --> K562
# treatment --> oxygen concentration 21%, 5%, 1%
oxygen_Dir = paste0(rawPath, "4_resources/public_CRISPR_raw_counts/K562_oxygen_jain/")
oxygen_raw_counts = vroom(paste0(oxygen_Dir, "raw_counts.tsv")) %>%
  as.matrix()
# load sample info 
oxygen_sampleinfo = vroom(paste0(oxygen_Dir, "sampleinfo.tsv"))
# more info
#readLines(paste0(oxygen_Dir, "README"))


### load sgRNAs library (Brunello)
grna_library = vroom(paste0(rawPath, "4_resources/crispr_libraries/brunello/original/broadgpp-brunello-library-contents_gRNAs.tsv")) %>%
  rename("gene" = "Target Gene Symbol") %>%
  select(sgRNA, gene, "sgRNA Target Sequence", "Target Context Sequence", "PAM Sequence", "Strand", "Genomic Sequence", "Position of Base After Cut (1-based)") %>%
  arrange(`Genomic Sequence`, `Position of Base After Cut (1-based)`)



### normalize count matrix 

## working only with APOBEC3A data for the moment

## varianceStabilizingTransformation
vst = varianceStabilizingTransformation(A3A_raw_counts)

# print normalized sample distributions
pdf(paste0(figDir, "normalization.pdf"), width=15, height=6)
par(mfrow=c(1,1))
boxplot(vst, outline=F, las=2, main="vst")
dev.off()


## create Biobase's 'ExpressionSets' (esets)
eset = new("ExpressionSet", exprs=vst)
pData(eset) = sampleinfo[match(colnames(eset), sampleinfo$sampleName),]
eset$group = factor(eset$group)
rownames(eset) = rownames(mat)
fData(eset) = grna_library[rownames(eset), ]

save(eset, file=paste0(dataDir, "eset.RData"))
#load(file=paste0(dataDir, "eset.RData"))


## PCAs
pr = prcomp(t(exprs(eset)))
pcvar = round(pr$sdev^2/sum(pr$sdev^2)*100, 2)

pdf(paste0(figDir, "pca.", i, ".pdf"), width=9, height=9)
par(mfrow=c(3, 3))
for(i in 2:3){
  plot(pr$x[,1], pr$x[,i], col=rainbow(5)[eset$group], xlab=paste0("PC1 (", pcvar[1], "%)"), ylab=paste0("PC", i, " (", pcvar[i], "%)"))
  text(pr$x[,1], pr$x[,i], labels=eset$sampleName, pos=4, cex=0.6, col=rainbow(5)[eset$group], srt=90)
}
plot(pr$x[,1], apply(exprs(eset), 2, median), col=rainbow(5)[eset$group], xlab=paste0("PC1 (", pcvar[1], "%)"), ylab="Col medians")
plot(pr$x[,1], colSums(mat), col=rainbow(5)[eset$group], xlab=paste0("PC1 (", pcvar[1], "%)"), ylab="Col sums")
plot(pr$x[,2], colSums(mat), col=rainbow(5)[eset$group], xlab=paste0("PC2 (", pcvar[2], "%)"), ylab="Col sums")
plot(pr$x[,3], colSums(mat), col=rainbow(5)[eset$group], xlab=paste0("PC3 (", pcvar[3], "%)"), ylab="Col sums")
plot(pr$x[,1], colMeans(mat == 0), col=rainbow(5)[eset$group], xlab=paste0("PC1 (", pcvar[1], "%)"), ylab="% zeros")
plot(pr$x[,2], colMeans(mat), col=rainbow(5)[eset$group], xlab=paste0("PC2 (", pcvar[2], "%)"), ylab="% zeros")
plot(pr$x[,3], colMeans(mat), col=rainbow(5)[eset$group], xlab=paste0("PC3 (", pcvar[3], "%)"), ylab="% zeros")
dev.off()



### differential lethality

sampleNames(eset) = eset$sampleName
eset$mouse[13:16] = 10:13
colnames(fData(eset))[2] = "gene"
eset = eset[!is.na(fData(eset)$gene),]
eset$treatment[13:14] = "NA"; eset$treatment[16] = "drug treatment"; eset$treatment = as.factor(eset$treatment)

# contrasts
conts =  mcp(group = c("tumor_untreated - tcbt_NA = 0",
                       "tumor_drugtreatment - tcbt_NA = 0",
                       "tumor_drugtreatment - tumor_untreated = 0"))

## run NB regression (loop across gene names)
allres = mclapply(unique(fData(eset)$gene), function(g){
  df = buildDF(g, eset, ref="tcbt_NA", cond="group")
  y = fitModel(df, conts, cond="group")
  m1 = tapply(df$counts, df$group, mean)
  list(df=df, res=y, mns=m1)
}, mc.cores=30)
names(allres) = unique(fData(eset)$gene)

save(allres, file=paste0(dataDir, "allres.RData"))

## parse results
cres = lapply(1:length(conts$group), function(i){
  y = do.call(rbind, lapply(allres, function(x) as.matrix(unname(x$res[i,, drop=F]))))
  colnames(y) = colnames(allres[[1]]$res)
  mns = do.call(rbind, lapply(allres, "[[", "mns"))
  y = cbind(y, mns)
  rownames(y) = names(allres)
  y
})
names(cres) = rownames(allres[[1]]$res)

cres = lapply(cres, function(x) {
  x = data.frame(x)
  x$pval.adj = p.adjust(x$pvals, method="BH")
  colnames(x) = sub("ci.", "", colnames(x))
  x
})

mns = cres[[1]][,5:9]

allres = cbind(mns, untrVsTCBT = cres[[1]][, c(1:4, 10)], treatVsTCBT = cres[[2]][, c(1:4, 10)], treatVsuntr = cres[[3]][, c(1:4, 10)])

# write tables
write.table(allres, file=paste0(tabDir, "results.xls"), quote=F, row.names=T, sep="\t")

for(i in names(cres)) write.table(cres[[i]], file=paste0(tabDir, i, ".results.xls"), quote=F, row.names=T, sep="\t")


## detect top 50 genes -sorted by max(|estimate|)- with significant regression (|estimate|>4 and FDR<0.05)
sel = unique(unlist(lapply(cres, function(x){
  x = x[order(abs(x$Estimate), decreasing=T),]
  x = x[x$pval.adj < 0.05 & abs(x$Estimate) > 4 ,]
  rownames(x)[1:min(nrow(x), 50)]
})))

mycols = c(1, '#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506')
mycols5 = mycols[c(1, 3, 5, 7, 9)]

system(paste0("mkdir ", figDir, "Stripcharts"))


## plots for the selected genes

# 2 barplots per gene
for(i in sel){
  
  y = buildDF(i, eset, ref="tcbt_NA")
  
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

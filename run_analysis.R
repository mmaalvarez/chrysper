
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
library(DESeq2)
library(multcomp)
library(MASS)
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
dataDir = paste0(mainDir, '/data/processed/') ; dir.create(dataDir, recursive = TRUE)

# output directories
outDir = paste0(mainDir, '/output/')
tabDir = paste0(outDir, 'tables/') ; dir.create(tabDir, recursive = TRUE)
figDir = paste0(outDir, 'figs/') ; dir.create(figDir, recursive = TRUE)

# temporary directory
tmpDir = "/local_scratch/malvarez/chrysper/tmp/"



#### load data

# raw data's root path
rawPath = "/g/strcombio/fsupek_data/CRISPR/"


### load tables of raw read counts and info from 3 different experiments

############################################
## our APOBEC3A experiment
# cell lines --> A549, LXF289, H358
# treatments --> DOX (A3A active) vs. control -- some also have ATRi
# altered genotypes --> TP53-KO, TP53-mut, HMCES-KO

A3A_Dir = paste0(rawPath, "2_processed_data/apobec/MAGeCK_MLE_analyses/raw_counts_combined/")

# load raw counts table
A3A_raw_counts = vroom(paste0(A3A_Dir, "A549_LXF289_H358__raw_counts.tsv")) %>%
  #sum counts across LANES of the same experiment
  group_by(sgRNA,Gene,sample_id,TP53,A3A,treatment,time) %>% 
  summarise_all(dplyr::funs(sum)) %>%
  
# load sample info 
A3A_sampleinfo = vroom(paste0(A3A_Dir, "sampleinfo.tsv"))

# more info
# see /g/strcombio/fsupek_data/CRISPR/1_raw_data/0_INFO/paths_info_A549_LXF289_H358.xlsx
############################################

############################################
## Travis lab's ATRi/TMZ experiment
# cell lines --> SW480, SW620, HT29
# treatments --> ATRi, TMZ vs. control
# altered genotypes --> ARID1A-KO, MSH6-KO

ATRi_TMZ_Dir = paste0(rawPath, "3_reports/stracker_lab/isabel_ATRi_TMZ/raw_counts_combined/")

# load raw counts table
ATRi_TMZ_raw_counts = vroom(paste0(ATRi_TMZ_Dir, "SW480_SW620_HT29__ARID1A_MSH6__raw_counts.tsv")) %>%
  #sum counts across LANES of the same experiment
  

# load sample info 
ATRi_TMZ_sampleinfo = vroom(paste0(ATRi_TMZ_Dir, "sampleinfo.tsv"))

# more info
#readLines("/g/strcombio/fsupek_data/CRISPR/1_raw_data/0_INFO/paths_info_isabel.txt")
############################################

############################################
## Jain et al oxygen experiment
# cell lines --> K562
# treatment --> oxygen concentration 21%, 5%, 1%

oxygen_Dir = paste0(rawPath, "4_resources/public_CRISPR_raw_counts/K562_oxygen_jain/")

# load raw counts table
oxygen_raw_counts = vroom(paste0(oxygen_Dir, "raw_counts.tsv"))

# load sample info 
oxygen_sampleinfo = vroom(paste0(oxygen_Dir, "sampleinfo.tsv"))

# more info
#readLines(paste0(oxygen_Dir, "README"))
############################################


### load sgRNAs library (Brunello)
grna_library = vroom(paste0(rawPath, "4_resources/crispr_libraries/brunello/original/broadgpp-brunello-library-contents_gRNAs.tsv")) %>%
  rename("gene" = "Target Gene Symbol") %>%
  # combine sgRNA and gene into an id
  unite(col = "id", sgRNA, gene, sep = "__", remove = F) %>%
  # keep interesting columns
  select(id, sgRNA, gene, "sgRNA Target Sequence", "Target Context Sequence", "PAM Sequence", "Strand", "Genomic Sequence", "Position of Base After Cut (1-based)") %>%
  # sort by chromosome and cut position
  arrange(`Genomic Sequence`, `Position of Base After Cut (1-based)`)



#### normalize count matrix 

### working only with APOBEC3A data for the moment

## varianceStabilizingTransformation --> normalizes read counts across samples
vst = A3A_raw_counts %>%
  # combine sgRNA and gene into an id
  unite(col = "id", sgRNA, gene, sep = "__") %>%
  # id to row names
  column_to_rownames("id") %>%
  # vst requires a matrix object
  as.matrix() %>%
  # run vst
  varianceStabilizingTransformation()


## print normalized sample distributions

# input format for ggplot
vst_plot_input =  vst %>%
  as.data.frame() %>%
  gather(key = sample,
         value = `vst-normalized sgRNA counts`) %>%
  # cell line / genotype variable
  mutate(group = sub("_.*", "", sample),
         group = sub("^[0-9]", "", group),
         group = sub("^[0-9]", "", group))
# fix sample order
vst_plot_input$sample = factor(vst_plot_input$sample,
                               ordered = TRUE,
                               levels = A3A_sampleinfo$sample)

## ggplot
boxplot = ggplot(data = vst_plot_input,
                 aes(x = sample,
                     y = `vst-normalized sgRNA counts`,
                     fill = factor(group))) +
  geom_boxplot(lwd = 0.1,
               outlier.size = 0.2,
               outlier.stroke = 0.2) +
  ggtitle("sgRNA read counts -- vst normalized -- A3A experiment") +
  xlab("") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 3,
                                   angle = 45,
                                   hjust = 1,
                                   vjust = 1.2),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(hjust = 0.5,
                                  size = 8),
        legend.position="none")

# save to file
ggsave(paste0(figDir,
              "vst_normalized_grna_counts_A3A.jpg"),
       plot = boxplot,
       device = "jpg",
       width = 10,
       height = 5.6,
       dpi = 300)



## store vst data as Biobase's 'ExpressionSets' (esets)
eset = new("ExpressionSet", exprs = vst)

# add sample info
pData(eset) = A3A_sampleinfo[match(colnames(eset), A3A_sampleinfo$sample),]
rownames(eset) = rownames(vst)
sampleNames(eset) = eset$sample

# set categorical variables as factors (unordered ('dummy' contrasts), except treatment)
eset$cell_line = factor(eset$cell_line)
eset$TP53 = factor(eset$TP53,
                   ordered = FALSE,
                   levels = c("wt", "mut", "KO"))
eset$HMCES = factor(eset$HMCES,
                    ordered = FALSE,
                    levels = c("wt", "KO"))
eset$A3A_vector = factor(eset$A3A_vector,
                         ordered = FALSE,
                         levels = c("no", "yes"))
eset$treatment = factor(eset$treatment,
                        ordered = TRUE,
                        levels = c("no", "DOX_IC25", "DOX_IC50", "DOX_ATRi"))
eset$replicate = factor(eset$replicate)

## set BACKWARD difference contrast matrices to use with ordered factors in the regression
# contrasts for treatment
contrast_matrix_four_levels =  matrix(c(-3/4,  1/4,  1/4, 1/4,
                                        -2/4, -2/4,  2/4, 2/4,
                                        -1/4, -1/4, -1/4, 3/4),
                                      ncol = 3,
                                      dimnames = list(c("no", "DOX_IC25", "DOX_IC50", "DOX_ATRi"),
                                                      c("DOX_IC25-vs-control", "DOX_IC50-vs-DOX_IC25", "DOX_ATRi-vs-DOX_IC50")))
contrasts(eset$treatment) <- contrast_matrix_four_levels


# add library info
fData(eset) = grna_library %>%
  filter(id %in% rownames(eset))


## save processed data
save(eset,
     file = paste0(dataDir, "eset_A3A.RData"))
# load it to save time
#load(file = paste0(dataDir, "eset_A3A.RData"))



### run NB regression (loop across gene names)
allres = mclapply(unique(fData(eset)$gene), function(g){
  df = buildDF(g, eset, ref="tcbt_NA", cond="group")
  y = fitModel(df, conts, cond="group")
  m1 = tapply(df$counts, df$group, mean)
  list(df=df, res=y, mns=m1)
}, mc.cores=30)
names(allres) = unique(fData(eset)$gene)

# save processed data
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

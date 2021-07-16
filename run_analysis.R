
#### load libraries
library(vroom)
library(magrittr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(DESeq2)
library(ggplot2)
library(gplots)
library(ggpubr)
library(MASS)
library(purrr)
library(broom)
library(stringr)
library(safejoin)
library(data.table)
library(readr)

# resolve function conflicts
library(conflicted)
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("melt", "data.table")
conflict_prefer("desc", "dplyr")



#### define and create directories

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




####### load dataset (tables of raw read counts and sample info), several options:


### a) testing dataset (few samples and genes)

## test raw counts table
raw_counts = vroom("test/test_raw_counts.tsv")

## test sample info 
sampleinfo = vroom("test/test_sampleinfo.tsv")

## test sgRNA library
grna_library = vroom("test/test_grna_library.tsv")



### b) 3 different experiments:

# raw data's root path
rawPath = "/g/strcombio/fsupek_data/CRISPR/"


######################################################################################################
## i) our APOBEC3A experiment
# cell lines --> A549, LXF289, H358
# treatments --> DOX (A3A active) vs. control -- some also have ATRi
# altered genotypes --> TP53-KO, TP53-mut, HMCES-KO

A3A_Dir = paste0(rawPath, "2_processed_data/apobec/MAGeCK_MLE_analyses/raw_counts_combined/")

# load raw counts table
raw_counts = vroom(paste0(A3A_Dir, "A549_LXF289_H358__raw_counts.tsv"))
  
# load sample info 
sampleinfo = vroom(paste0(A3A_Dir, "sampleinfo.tsv"))

# more info
# see /g/strcombio/fsupek_data/CRISPR/1_raw_data/0_INFO/paths_info_A549_LXF289_H358.xlsx
######################################################################################################


######################################################################################################
## ii) Travis lab's ATRi/TMZ experiment
# cell lines --> SW480, SW620, HT29
# treatments --> ATRi, TMZ vs. control
# altered genotypes --> ARID1A-KO, MSH6-KO

ATRi_TMZ_Dir = paste0(rawPath, "3_reports/stracker_lab/isabel_ATRi_TMZ/raw_counts_combined/")

# load raw counts table
raw_counts = vroom(paste0(ATRi_TMZ_Dir, "SW480_SW620_HT29__ARID1A_MSH6__raw_counts.tsv"))

# load sample info 
sampleinfo = vroom(paste0(ATRi_TMZ_Dir, "sampleinfo.tsv"))

# more info
#readLines("/g/strcombio/fsupek_data/CRISPR/1_raw_data/0_INFO/paths_info_isabel.txt")
######################################################################################################


######################################################################################################
## iii) Jain et al oxygen experiment
# cell lines --> K562
# treatment --> oxygen concentration 21%, 5%, 1%

oxygen_Dir = paste0(rawPath, "4_resources/public_CRISPR_data/1_raw_data/K562_oxygen_jain/")

# load raw counts table
raw_counts = vroom(paste0(oxygen_Dir, "raw_counts.tsv"))

# load sample info 
sampleinfo = vroom(paste0(oxygen_Dir, "sampleinfo.tsv"))

# more info
#readLines(paste0(oxygen_Dir, "README"))
######################################################################################################


# ### (all, NOT USED) load gene essentiality (mean DEMETER2 score <0 --> gene inhibition decreases fitness consistently across cell lines)
# D2_score = vroom(paste0(rawPath, "4_resources/genes_info/essentiality/depmap/demeter/mean_D2.tsv"))

### (all) load sgRNAs library (Brunello)
grna_library = vroom(paste0(rawPath, "4_resources/crispr_libraries/brunello/original/broadgpp-brunello-library-contents_gRNAs.tsv")) %>%
  rename("gene" = "Target Gene Symbol") %>%
  # update Non-Targeting sgRNAs name (remove trailing " Control" part)
  mutate(gene = gsub(" Control$", "", gene)) %>%
  # combine sgRNA and gene into an id
  unite(col = "sgRNA_id", sgRNA, gene, sep = "__", remove = F) %>%
  # keep interesting columns
  select(sgRNA_id, sgRNA, gene, "sgRNA Target Sequence", "Target Context Sequence", "PAM Sequence", "Strand", "Genomic Sequence", "Position of Base After Cut (1-based)") %>%
  # sort by chromosome and cut position
  arrange(`Genomic Sequence`, `Position of Base After Cut (1-based)`)
  # add demeter2 score
  #%>% merge(D2_score)



#### normalize count matrix (to allow comparisons between samples)

### working only with APOBEC3A data for the moment

## option a): varianceStabilizingTransformation --> normalizes read counts across samples
vst = raw_counts %>%
  # combine sgRNA and gene into an id
  unite(col = "sgRNA_id", sgRNA, gene, sep = "__") %>%
  # sgRNA_id to row names
  column_to_rownames("sgRNA_id") %>%
  # vst requires a matrix object
  as.matrix() %>%
  # run vst
  varianceStabilizingTransformation()

## print normalized sample distributions

# input format required for ggplot
vst_plot_input =  vst %>%
  as.data.frame() %>%
  gather(key = sample,
         value = `vst-normalized sgRNA counts`) %>%
  # cell line / genotype variable
  mutate(group = sub("_.*", "", sample),
         group = sub("^[0-9]", "", group),
         group = sub("^[0-9]", "", group))
# set a fixed sample order
vst_plot_input$sample = factor(vst_plot_input$sample,
                               ordered = TRUE,
                               levels = sampleinfo$sample)
## boxplot using ggplot
boxplot = ggplot(data = vst_plot_input,
                 aes(x = sample,
                     y = `vst-normalized sgRNA counts`,
                     fill = factor(group))) +
  geom_boxplot(lwd = 0.1,
               outlier.size = 0.2,
               outlier.stroke = 0.2) +
  #ggtitle("[A3A, ATRi/TMZ, oxygen] experiment samples") +
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
ggsave(paste0(figDir, "vst_normalized_grna_counts.jpg"),
       plot = boxplot, device = "jpg", width = 10, height = 5.6, dpi = 300)


## option b) (preferred here): as Brunello does contain control non-targeting sgRNAs, we can simply use the ln(sum(Non-Targeting sgRNA counts)) as an offset in the NB regression
# store this offset in a table for later adding it to the data frames created in buildDF()
raw_counts_offset = raw_counts %>%
  # long format
  gather(key = "sample",
         value = "counts",
         -c(sgRNA, gene)) %>%
  group_by(sample) %>%
  # create offset column
  mutate(ln_sum_nontargeting = log(sum(ifelse(gene == "Non-Targeting",
                                              yes = counts,
                                              no = 0),
                                       na.rm = TRUE))) %>%
  select(sample, ln_sum_nontargeting) %>%
  distinct()

# prepare raw counts table for conversion into eset
raw_counts_for_eset = raw_counts %>%
  # combine sgRNA and gene into an id
  unite(col = "sgRNA_id", sgRNA, gene, sep = "__") %>%
  # sgRNA_id to row names
  column_to_rownames("sgRNA_id") %>%
  # ExpressionSet requires a matrix object as input
  as.matrix()


## Other normalization alternatives:
# 1- MAGeCK median-normalization --> yields table of normalized counts
# 2- If the library contains control sgRNAs (non-targeting (like in Brunello), or safe-targeting)
#    i) MAGeCK median-normalization using only this control set as reference
# 3- rlog() (another 'DESeq2' function)




#### store data as Biobase's 'ExpressionSets' (esets): a Class to Contain and Describe High-Throughput Expression Level Assays

eset = new("ExpressionSet", exprs = raw_counts_for_eset) # in case of wanting to use the vst normalized data, use 'exprs = vst'

# which is the order of sample names in eset
sample_order = colnames(eset)

# which is the order of sgRNA ids in eset
grna_order = rownames(eset)


### add sample info to eset's 'phenoData' section

# make sure order of sampleinfo's rows is the same as in 'sample_order'
# only samples that are BOTH in the counts AND sampleinfo tables are included (in theory, all of them)
sampleinfo_sorted = left_join(data.frame(sample = sample_order),
                              sampleinfo,
                              by = "sample")

pData(eset) = sampleinfo_sorted


## set the phenoData's conditions that are categorical variables as factors: either unordered (use 'dummy' contrasts by default), or unordered (treatment and time, set contrasts below)
eset$cell_line = factor(eset$cell_line,
                        ordered = FALSE,
                        levels = c("A549", "LXF289", "H358"))
eset$TP53 = factor(eset$TP53,
                   ordered = FALSE,
                   levels = c("wt", "mut", "KO"))
eset$HMCES = factor(eset$HMCES,
                    ordered = FALSE,
                    levels = c("wt", "KO"))
eset$A3A_vector = factor(eset$A3A_vector,
                         ordered = FALSE,
                         levels = c("no", "yes"))
eset$replicate = factor(eset$replicate,
                        ordered = FALSE,
                        levels = c(1, 2))
# recategorize treatment variable (treatment_recat)
pData(eset) %<>%
  mutate(treatment_recat = ifelse(treatment == "no",
                                  yes = "control",
                                  no = ifelse(treatment == "DOX_IC25",
                                              yes = "DOX_IC25",
                                              no = ifelse(treatment %in% c("DOX_IC50", "DOX_ATRi"),
                                                          yes = "DOX_IC50orATRi",
                                                          no = "Warning: treatment other than no,DOX_IC25,DOX_IC50,DOX_ATRi"))),
         # set treatment_recat as ordered factor, and levels order
         treatment_recat = factor(treatment_recat,
                                  ordered = TRUE,
                                  levels = c("control",
                                             "DOX_IC25",
                                             "DOX_IC50orATRi")))
# categorize time variable (time_cat)
pData(eset) %<>%
  mutate(time_cat = ifelse(time == 0,
                           yes = "time_0",
                           no = ifelse(time %in% c(3,5,6,9,10,12,13),
                                       yes = "time_3_5_6_9_10_12_13",
                                       no = ifelse(time %in% c(15,17),
                                                   yes = "time_15_17",
                                                   no = "Warning: time point other than 0,3,5,6,9,10,12,13,15,17"))),
         # set time_cat as ordered factor, and levels order
         time_cat = factor(time_cat,
                           ordered = TRUE,
                           levels = c("time_0",
                                      "time_3_5_6_9_10_12_13",
                                      "time_15_17")))

## set BACKWARD difference contrast matrices to use with ordered factors in the regression
# contrasts for treatment_recat
contrast_matrix_three_levels_treatment_recat = matrix(c(-2/3, 1/3, 1/3,
                                                        -1/3,-1/3, 2/3),
                                                      ncol = 2,
                                                      dimnames = list(c("control",
                                                                        "DOX_IC25",
                                                                        "DOX_IC50orATRi"),
                                                                      c("DOX_IC25-vs-control", "DOX_IC50orATRi-vs-DOX_IC25")))
contrasts(eset$treatment_recat) <- contrast_matrix_three_levels_treatment_recat

# contrasts for time_cat
contrast_matrix_three_levels_time_cat = matrix(c(-2/3, 1/3, 1/3,
                                                 -1/3,-1/3, 2/3),
                                               ncol = 2,
                                               dimnames = list(c("time_0",
                                                                 "time_3_5_6_9_10_12_13",
                                                                 "time_15_17"),
                                                               c("middle-vs-zero", "late-vs-middle")))
contrasts(eset$time_cat) <- contrast_matrix_three_levels_time_cat

## list of conditions (column names in 'sampleinfo'):
# "cell_line","TP53","HMCES","A3A_vector","treatment","time","replicate","treatment_recat","time_cat"
conditions = names(pData(eset))[! names(pData(eset)) %in% c("sample", "id")]


### add library info to eset's 'featureData' section

# make sure order of grna_library's rows is the same as in 'grna_order'
# only sgRNAs that are BOTH in the counts AND library tables are included (in theory, all of them)
grna_library_sorted = left_join(data.frame(sgRNA_id = grna_order),
                                grna_library,
                                by = "sgRNA_id")

fData(eset) = grna_library_sorted


## phenoData's 'sampleNames' section are just numbers -- replace with the actual sample names (SHOULD BE IN THE SAME ORDER)
sampleNames(eset) = eset$sample


## save processed data
save(eset,
     file = paste0(dataDir, "eset.RData"))
# load it to save time
#load(file = paste0(dataDir, "eset.RData")) ; conditions = names(pData(eset))[! names(pData(eset)) %in% c("sample", "id")]




#### NB regression

### load NB regression custom functions
source("src/functions.R")

### run NB regression (loop across genes) --> for full A3A dataset, it takes ~30 min. in fsupeksvr
gene_names = fData(eset) %>%
  select(gene) %>%
  filter(gene != "Non-Targeting") %>%
  distinct() %>%
  pull(gene)

NBres = mclapply(gene_names, function(g){
  
  # build dataframe (for a given gene) for regression
  df = buildDF(g, eset, conditions, raw_counts_offset)
  
  # run NB regression on that gene
  y = fitModel(g, eset, conditions, df)
  
}, mc.cores = 30)

# name elements by their analyzed gene's name
names(NBres) = gene_names

## save processed data
save(NBres, file = paste0(dataDir, "NBres.RData"))
#load(file = paste0(dataDir, "NBres.RData"))




#### parse NB results

## merge all genes results (i.e. the elements from NBres list) into megatable 'allres'

# apply broom::tidy() to each gene from NBres
NBres_tidy = NBres %>%
  map(., ~tidy(.x$model)) %>%
  # remove sgRNA rows
  map(., ~filter(.x, ! str_detect(term, "sgRNAs")))

# now merge NB results for all genes
allres = eat(NBres_tidy[[1]],
             NBres_tidy[-1],
             .by = "term") %>%
  # the first dataset's variables remained unchanged, so append dataset name to them
  rename_with(.fn = ~paste0(names(NBres_tidy)[1],
                            "_estimate"),
              .cols = all_of("estimate")) %>%
  rename_with(.fn = ~paste0(names(NBres_tidy)[1],
                            "_std.error"),
              .cols = all_of("std.error")) %>%
  rename_with(.fn = ~paste0(names(NBres_tidy)[1],
                            "_statistic"),
              .cols = all_of("statistic")) %>%
  rename_with(.fn = ~paste0(names(NBres_tidy)[1],
                            "_p.value"),
              .cols = all_of("p.value")) %>%
  # long format, by term
  gather(key = "estimate.SE.stat.pval",
         value = "value",
         -term) %>%
  # split 'estimate.SE.stat.pval' into 'estimate.SE.stat.pval' and 'gene' columns
  separate(estimate.SE.stat.pval,
           c("gene",
             "estimate.SE.stat.pval"),
           sep = "\\_") %>%
  # one column estimate, one column SE, another statistic, another p.value
  pivot_wider(names_from = c(estimate.SE.stat.pval),
                     values_from = c(value)) %>%
  # remove NA rows
  na.omit() %>%
  # remove intercepts
  filter(term != "(Intercept)") %>%
  # add column with the condition tested (e.g. treatment_recatDOX_IC25-vs-control.. is "treatment") OR if it is interaction (contains ":")
  mutate(condition = ifelse(str_detect(term, "\\:"),
                            yes = "interaction",
                            no = gsub("_.*",
                                       "",
                                       term))) %>%
  # correct p-values (FDR)
  mutate(fdr = p.adjust(p.value,
                        method = "BH")) %>%
  # order and sort
  select(gene, condition, term, estimate, std.error, statistic, p.value, fdr) %>%
  arrange(fdr, p.value)
  
## write table to file
write_tsv(allres,
          paste0(tabDir, "results.tsv"))



#### detect top gene hits

hits = allres %>%
  # some filters for significance, e.g. |estimate|>4 & FDRs<0.05
  filter(abs(estimate) > 0.001 &
           fdr < 0.95) %>%
  # sorted by e.g. descending |estimate|
  arrange(desc(abs(estimate)))


## plot hits

for(gene in hits){
  
  # plot sgRNA real (raw) counts per gene×treatment×time
  p_real_counts = ggplot(data = NBres[[gene]]$df,
                         aes(x = treatment_recat,
                             y = `real counts`,
                             colour = sgRNA)) +
                    geom_boxplot() +
                    facet_grid(facets = ~time_cat) +
                    theme_bw() +
                    ggtitle(gene)
  ggsave(paste0(figDir, gene, "_raw_grna_counts.jpg"),
         plot = p_real_counts, device = "jpg", width = 10, height = 5.6, dpi = 300)

  # plot real (raw) vs. fitted counts per gene×experiment_id, we want them to be similar (stratify by experiment group "id")
  p_real_vs_fitted_counts = ggplot(data = NBres[[gene]]$df,
                                   aes(x = `real counts`,
                                       y = `fitted counts`,
                                       group = id)) +
                              geom_point() +
                              geom_smooth(method = lm) + 
                              stat_cor(method = "pearson",
                                       alternative = "two.sided",
                                       cor.coef.name = "R") + 
                              facet_grid(facets = ~id) + 
                              ylab("NB fitted counts") + 
                              theme_minimal() +
                              ggtitle(gene)
  ggsave(paste0(figDir, gene, "_raw_vs_fitted_grna_counts.jpg"),
         plot = p_real_vs_fitted_counts, device = "jpg", width = 10, height = 5.6, dpi = 300)
}

# heatmap sample×sgRNAs
pdf(paste0(figDir, "top_hits_heatmap.pdf"), width=10, height=5.6)
heatmap.2(exprs(eset[fData(eset)$gene %in% hits, ]),
          trace="none", scale="row",
          ColSideColors=c("gray", "yellow", "red")[eset$treatment_recat])
dev.off()

gene = "A1BG"
NBres[[gene]]$df %>%
  mutate(treatment_2_levels = ifelse(treatment == "no",
                                     yes = "control",
                                     no = "treated")) %>%
  ggplot(aes(x = time_cat,
           y = `fitted counts`,
           color = treatment_2_levels)) +
  geom_violin(width = 0.2,
              position = position_dodge(width = 0.3)) +
  geom_boxplot(width = 0.1,
               position = position_dodge(width = 0.3)) +
  stat_summary(aes(group = treatment_2_levels),
               geom = "point",
               fun = "mean",
               position = position_dodge(width = 0.3)) +
  stat_summary(aes(group = treatment_2_levels),
               geom = "line",
               fun = "mean",
               position = position_dodge(width = 0.3)) +
  ggtitle(paste0("gene: ", gene)) +
  theme_minimal()
  



data_glm = NBres[["A1BG"]]$df
# plot number of counts vs time including the regression
p_counts_vs_time = ggplot(data_glm[data_glm$treatment_recat == 'control', ],
                          aes(x = time_cat,
                              y = `real counts`)) + 
  geom_boxplot()



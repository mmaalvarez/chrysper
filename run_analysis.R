#### load libraries
library(vroom)
library(magrittr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(Biobase)
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
library(effects)          # plot interaction terms
library(readr)
library(safejoin)
library(devtools)
library(sjmisc)           # for plot_model function
library(sjPlot)           # for plot_model function
library(ggiraphExtra)
library(insight)
library(moonBook)
library(ggeffects)
library(Biobase)


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
raw_counts_t = vroom("test/test_raw_counts.tsv")

## test sample info 
sampleinfo_t = vroom("test/test_sampleinfo.tsv")

## test sgRNA library
grna_library_t = vroom("test/test_grna_library.tsv")



### b) 3 different experiments:

# raw data's root path
rawPath = "/g/strcombio/fsupek_data/CRISPR/"


######################################################################################################
## i) our APOBEC3A experiment (A3A)
# cell lines --> A549, LXF289, H358
# treatments --> DOX (A3A active) vs. control -- some also have ATRi
# altered genotypes --> TP53-KO, TP53-mut, HMCES-KO

A3A_Dir = paste0(rawPath, "2_processed_data/apobec/MAGeCK_MLE_analyses/raw_counts_combined/")

# load raw counts table
raw_counts_A3A = vroom(paste0(A3A_Dir, "A549_LXF289_H358__raw_counts.tsv"))

# load sample info 
sampleinfo_A3A = vroom(paste0(A3A_Dir, "sampleinfo.tsv"))

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


oxygen_Dir = "/g/strcombio/fsupek_cancer3/malvarez/public_CRISPR_data/1_raw_data/K562_oxygen_jain/"

# load raw counts table
raw_counts = vroom(paste0(oxygen_Dir, "raw_counts.tsv")) # overwriting variables?

# load sample info 
sampleinfo = vroom(paste0(oxygen_Dir, "sampleinfo.tsv"))

# more info
#readLines(paste0(oxygen_Dir, "README"))
######################################################################################################


# ### (all, NOT USED) load gene essentiality (mean DEMETER2 score <0 --> gene inhibition decreases fitness consistently across cell lines)
# D2_score = vroom(paste0(rawPath, "4_resources/genes_info/essentiality/depmap/demeter/mean_D2.tsv"))

### (all) load sgRNAs library (Brunello)
grna_library = vroom(paste0(rawPath, "4_resources/crispr_libraries/brunello/original/broadgpp-brunello-library-contents_gRNAs.tsv")) %>%
  # update Non-Targeting sgRNAs name (remove trailing " Control" part)
  mutate(gene = gsub(" Control$","", gene)) %>%
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

## as Brunello does contain control non-targeting sgRNAs, we can simply use the ln(sum(Non-Targeting sgRNA counts)) as an offset in the NB regression
# store this offset in a table for later adding it to the data frames created in buildDF()
raw_counts_offset = raw_counts_A3A %>%
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
raw_counts_for_eset = raw_counts_A3A %>%
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
# 3- rlog() or varianceStabilizingTransformation() ('DESeq2' functions)




#### store data as Biobase's 'ExpressionSets' (esets): a Class to Contain and Describe High-Throughput Expression Level Assays

eset = new("ExpressionSet", exprs = raw_counts_for_eset) # in case of wanting to use the vst normalized data, use 'exprs = vst'

# which is the order of sample names in eset
sample_order = colnames(eset)

# which is the order of sgRNA ids in eset
grna_order = rownames(eset) # takes the first column


### add sample info to eset's 'phenoData' section

# make sure order of sampleinfo's rows is the same as in 'sample_order'
# only samples that are BOTH in the counts AND sampleinfo tables are included (in theory, all of them)
sampleinfo_sorted = left_join(data.frame(sample = sample_order),
                              sampleinfo_A3A,
                              by = "sample")

pData(eset) = sampleinfo_sorted # columns: sample, id, cell_line, TP53, HMCES, A3A_vector, treatment, time, replicate


## set the phenoData's conditions that are categorical variables as factors: either unordered (use 'dummy' contrasts by default), or unordered (treatment and time, set contrasts below)
eset$cell_line = factor(eset$cell_line,
                        ordered = FALSE,
                        levels = c("A549", "LXF289", "H358"))     #for the first experiment
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

## recategorize treatment variable into ordered factor (treatment_recat)
# a) TEST A3A dataset
pData(eset) %<>%
  mutate(treatment_recat = ifelse(treatment == "no",
                                  yes = "control",
                                  no = "treated"),
         # set treatment_recat as ordered factor, and levels order
         treatment_recat = factor(treatment_recat,
                                  ordered = TRUE,
                                  levels = c("control",
                                             "treated")))
# b) FULL A3A DATASET --> treated only if A3A + DOX 
pData(eset) %<>%
  mutate(treatment_recat = ifelse((A3A_vector == 'no' & treatment == 'no') |
                                    (A3A_vector == 'yes' & treatment == 'no') |
                                    (A3A_vector == 'no' & treatment == 'DOX_IC25') |
                                    (A3A_vector == 'no' & treatment == 'DOX_ATRi'),
                                  yes = "control",
                                  no = "treated"),
         # set treatment_recat as ordered factor, and levels order
         treatment_recat = factor(treatment_recat,
                                  ordered = TRUE,
                                  levels = c("control",
                                             "treated")))

# categorize time variable into ordered factor (time_cat)
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

## set BACKWARD difference contrast matrices to use with time_cat in the regression

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

# IMPORTANT: depending on the analysis, change formula in src/functions.R

### load NB regression custom functions
source("src/functions.R")

### run NB regression (loop across genes) --> for full A3A dataset, it takes 8? hours in fsupeksvr

gene_names = fData(eset) %>%
  select(gene) %>%
  filter(gene != "Non-Targeting") %>%
  distinct() %>%
  pull(gene)

NBres = mclapply(gene_names, function(g){
  
  # build dataframe (for a given gene) for regression
  df = buildDF(g, eset, conditions, raw_counts_offset)
  
  ## this whole step should probably be moved to the eset preparation step
  df = df %>%
    # add time_0 counts to treated samples (i.e. t0 is shared for both control and treated samples) 
    filter(time_cat == 'time_0') %>%
    mutate(treatment_recat = replace(treatment_recat, treatment_recat == 'control','treated')) %>%
    rbind(df)
  # redefine contrast for time_cat to add the new time_0 treated rows
  contrasts(df$time_cat) <- contrast_matrix_three_levels_time_cat
  ##
  
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
  map(., ~tidy(.x$model))

# now merge NB results for all genes into a megatable (takes time)
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

# to detect genes with just negative slope for treated samples
hits = data.frame(matrix(nrow=0,
                         ncol=length(colnames(allres))))
names(hits) = colnames(allres)
hits = as_tibble(hits)
for (g in unique(allres$gene)){
  treatment_L = allres %>% filter(gene == g & term == 'treatment_recat.L')
  time_mid_0 = allres %>% filter(gene == g & term == 'time_catmiddle-vs-zero')
  time_late_mid = allres %>% filter(gene == g & term == 'time_catlate-vs-middle')
  interaction_mid_0 = allres %>% filter(gene == g & term == 'time_catmiddle-vs-zero:treatment_recat.L')
  interaction_late_mid = allres %>% filter(gene == g & term == 'time_catlate-vs-middle:treatment_recat.L')
  ifelse(interaction_mid_0$fdr < 0.05 ||
           interaction_late_mid$fdr < 0.05,
         yes = (hits = rbind(hits, interaction_mid_0,
                             interaction_late_mid,
                             treatment_L,
                             time_late_mid,
                             time_mid_0)),
         no = (hits = hits))
}


## plot hits

for(gene in unique(hits$gene)){
  data_df = NBres[[gene]]$df
  data_model = NBres[[gene]]$model
  new_model_names = c('counts', 'time_cat', 'treatment_recat', 'ln_sum_nontargeting')
  names(data_model$model) = new_model_names
  p_model = plot_model(data_model, 
                       type = "pred", 
                       terms = c('time_cat', 'treatment_recat'),
                       value.offset = 'ln_sum_nontargeting') +
    geom_line() +
    ggtitle(gene) +
    theme_minimal()
  show(p_model)
  ggsave(path = figDir, filename = paste0(gene, "_counts_vs_time.jpg"),
          plot = p_model, device = "jpg", width = 10, height = 5.6, dpi = 300)  
}

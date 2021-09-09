#### load libraries
library(vroom)
library(magrittr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(DESeq2)
library(ggplot2)
# library(gplots)
# library(ggpubr)
library(MASS)
library(purrr)
# library(broom)
# library(stringr)
# library(safejoin)
library(data.table)
# library(effects)          # plot interaction terms
# library(readr)
# library(safejoin)
# library(devtools)
# library(sjmisc)           # for plot_model function
# library(sjPlot)           # for plot_model function
# library(ggiraphExtra)
# library(insight)
# library(moonBook)
# library(ggeffects)
# library(Biobase)
# 

# resolve function conflicts
library(conflicted)
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
#conflict_prefer("melt", "data.table")
#conflict_prefer("desc", "dplyr")



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



####### load datasets (tables of raw read counts and sample info), several options:

# raw data's root path
rawPath = "/g/strcombio/fsupek_data/CRISPR/"

### load sgRNAs library (Brunello), common for the 3 datasets
grna_library = vroom(paste0(rawPath, "4_resources/crispr_libraries/brunello/original/broadgpp-brunello-library-contents_gRNAs.tsv")) %>%
  # update Non-Targeting sgRNAs name (remove trailing " Control" part)
  mutate(gene = gsub(" Control$","", gene)) %>%
  # combine sgRNA and gene into an id
  unite(col = "sgRNA_id", sgRNA, gene, sep = "__", remove = F) %>%
  # keep interesting columns
  select(sgRNA_id, sgRNA, gene, "sgRNA Target Sequence", "Target Context Sequence", "PAM Sequence", "Strand", "Genomic Sequence", "Position of Base After Cut (1-based)") %>%
  # sort by chromosome and cut position
  arrange(`Genomic Sequence`, `Position of Base After Cut (1-based)`)

### raw data from 3 different experiments:

######################################################################################################
## i) our APOBEC3A experiment (A3A)
# cell lines --> A549, LXF289, H358
# treatments --> DOX (A3A active) vs. control -- some also have ATRi
# BUT H358 does not have control samples
# altered genotypes --> TP53-KO, TP53-mut, HMCES-KO
# replicates --> A549 TP53-/- A3A+DOX, but it is weird, maybe not include

A3A_Dir = paste0(rawPath, "2_processed_data/apobec/MAGeCK_MLE_analyses/raw_counts_combined/")

# load raw counts table
raw_counts_A3A = vroom(paste0(A3A_Dir, "A549_LXF289_H358__raw_counts.tsv"))

# load sample info 
sampleinfo_A3A = vroom(paste0(A3A_Dir, "sampleinfo.tsv"))

# more info --> /g/strcombio/fsupek_data/CRISPR/1_raw_data/0_INFO/paths_info_A549_LXF289_H358.xlsx
######################################################################################################


######################################################################################################
## ii) Travis lab's ATRi/TMZ experiment
# cell lines --> SW480, SW620, HT29
# treatments --> ATRi, TMZ vs. control
# altered genotypes --> ARID1A-KO, MSH6-KO

ATRi_TMZ_Dir = paste0(rawPath, "3_reports/stracker_lab/isabel_ATRi_TMZ/raw_counts_combined/")

# load raw counts table
raw_counts_ATRi_TMZ = vroom(paste0(ATRi_TMZ_Dir, "SW480_SW620_HT29__ARID1A_MSH6__raw_counts.tsv"))

# load sample info 
sampleinfo_ATRi_TMZ = vroom(paste0(ATRi_TMZ_Dir, "sampleinfo.tsv"))

# more info --> readLines("/g/strcombio/fsupek_data/CRISPR/1_raw_data/0_INFO/other_labs/paths_info_stracker.txt")
######################################################################################################


######################################################################################################
## iii) Jain et al oxygen experiment
# cell lines --> K562
# treatment --> oxygen concentration 21%, 5%, 1%

oxygen_Dir = "/g/strcombio/fsupek_cancer3/malvarez/public_CRISPR_data/1_raw_data/K562_oxygen_jain/"

# table of genes whose name differs in the other datasets
gene_names_update = vroom(paste0(oxygen_Dir, "genes_different_names_other_datasets.tsv")) %>%
  rename("gene" = "gene_name_oxygen_dataset")

# load raw counts table
raw_counts_oxygen = vroom(paste0(oxygen_Dir, "raw_counts.tsv")) %>%
  ## compatibility with the other datasets
  # change some gene names
  mutate(gene = gsub("NO_CURRENT_.*", "Non-Targeting", gene)) %>%
  merge(gene_names_update, all = T) %>%
  mutate(gene = ifelse(is.na(gene_name_other_datasets),
                       gene,
                       gene_name_other_datasets)) %>%
  select(-gene_name_other_datasets) %>%
  # convert sgRNA sequences into sgRNA ids
  rename("sgRNA Target Sequence" = "sgRNA") %>%
  merge(select(grna_library, sgRNA, `sgRNA Target Sequence`),
        by = "sgRNA Target Sequence") %>%
  select(-`sgRNA Target Sequence`) %>%
  relocate(sgRNA, .before = "gene")

# load sample info 
sampleinfo_oxygen = vroom(paste0(oxygen_Dir, "sampleinfo.tsv"))

# more info --> readLines(paste0(oxygen_Dir, "README"))
######################################################################################################


### combine datasets
raw_counts_3_datasets = merge(raw_counts_A3A, raw_counts_ATRi_TMZ) %>%
  merge(raw_counts_oxygen)

sampleinfo_3_datasets = merge(sampleinfo_A3A, sampleinfo_ATRi_TMZ, all = T) %>%
  merge(sampleinfo_oxygen, all = T) %>%
  # control or treated
  mutate(type = ifelse(! is.na(A3A_vector), 
                       # if it is a dataset with A3A, we'll consider a sample to be control if noA3A and noATRi -- A3A without DOX is weird, and noA3A with ATRi is treated...
                       ifelse((A3A_vector == 'no' & treatment == 'no') |         # "true control"
                                (A3A_vector == 'no' & treatment == 'DOX_IC25') | # "pseudo-replicate 1"
                                (A3A_vector == 'no' & treatment == 'DOX_IC50'),  # "pseudo-replicate 2"
                              "control",
                              "treated"),
                       # other datasets
                       ifelse(treatment == "no",
                              "control",
                              "treated"))) %>%
  # also consider each special genotype as an independent cell line
  mutate(cell_line = ifelse(! is.na(TP53),
                            paste0(cell_line, "_TP53", TP53),
                            cell_line),
         cell_line = ifelse(! is.na(HMCES),
                            paste0(cell_line, "_HMCES", HMCES),
                            cell_line),
         cell_line = ifelse(! is.na(ARID1A),
                            paste0(cell_line, "_ARID1A", ARID1A),
                            cell_line),
         cell_line = ifelse(! is.na(MSH6),
                            paste0(cell_line, "_MSH6", MSH6),
                            cell_line)) %>%
  select(sample, cell_line, time, type, replicate) %>%
  arrange(cell_line, time, type, replicate)

sampleinfo_3_datasets %>%
  filter(type == "control") %>%
  select(cell_line, time) %>%
  table()
#                   time 0 3 5 6 7 9 10 12 14 15
#-----------------------------------------------
# A549_TP53KO_HMCESwt    1 2 0 2 0 2  0  2  0  2
# A549_TP53wt_HMCESwt    1 2 0 2 0 2  0  2  0  2
# HT29_ARID1Awt_MSH6KO   1 0 0 0 1 0  0  0  1  0
# HT29_ARID1Awt_MSH6wt   1 0 0 0 1 0  0  0  1  0
# K562                   3 0 0 0 0 2  0  0  0  3
# LXF289_TP53mut_HMCESwt 1 0 3 0 0 0  3  0  0  3
# SW480_ARID1AKO_MSH6wt  1 0 0 0 1 0  0  0  1  0
# SW480_ARID1Awt_MSH6wt  1 0 0 0 1 0  0  0  1  0
# SW620_ARID1AKO_MSH6wt  2 0 0 0 2 0  0  0  2  0
# SW620_ARID1Awt_MSH6KO  1 0 0 0 1 0  0  0  1  0
# SW620_ARID1Awt_MSH6wt  3 0 0 0 3 0  0  0  3  0


#### normalize count matrix (to allow comparisons between samples)

## as Brunello contains control non-targeting sgRNAs, we can simply use the ln(sum(Non-Targeting sgRNA counts)) as an offset in the NB regression
# store this offset in a table for later adding it to the data frames created in buildDF()
raw_counts_offset = raw_counts_3_datasets %>%
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
raw_counts_for_eset = raw_counts_3_datasets %>%
  # combine sgRNA and gene into an id
  unite(col = "sgRNA_id", sgRNA, gene, sep = "__") %>%
  # sgRNA_id to row names
  column_to_rownames("sgRNA_id") %>%
  # ExpressionSet requires a matrix object as input
  as.matrix()


#### store data as Biobase's 'ExpressionSets' (esets): a Class to Contain and Describe High-Throughput Expression Level Assays

eset = new("ExpressionSet", exprs = raw_counts_for_eset) # in case of wanting to use the vst normalized data, use 'exprs = vst'

# which is the order of sample names in eset
sample_order = colnames(eset)

# which is the order of sgRNA ids in eset
grna_order = rownames(eset) # takes the first column


### add sample info to eset's 'phenoData' section

# make sure order of sampleinfo's rows is the same as in the eset, i.e. 'sample_order'
# only samples that are BOTH in the counts AND sampleinfo tables are included (in theory, all of them)
sampleinfo_sorted = left_join(data.frame(sample = sample_order),
                              sampleinfo_3_datasets,
                              by = "sample")

pData(eset) = sampleinfo_sorted # columns: sample, id, cell_line, TP53, HMCES, A3A_vector, treatment, time, replicate


## set the phenoData's conditions that are categorical variables as factors: either unordered (use 'dummy' contrasts by default), or unordered (treatment and time, set contrasts below)
eset$cell_line = factor(eset$cell_line,
                        ordered = FALSE,
                        levels = c("A549_TP53KO_HMCESwt",
                                   "A549_TP53wt_HMCESwt",
                                   "LXF289_TP53mut_HMCESKO",
                                   "LXF289_TP53mut_HMCESwt",
                                   "H358_TP53mut_HMCESwt",   # this cell line does not have control samples
                                   "HT29_ARID1Awt_MSH6KO",
                                   "HT29_ARID1Awt_MSH6wt",
                                   "SW480_ARID1AKO_MSH6wt",
                                   "SW480_ARID1Awt_MSH6wt",
                                   "SW620_ARID1AKO_MSH6wt",
                                   "SW620_ARID1Awt_MSH6KO",
                                   "SW620_ARID1Awt_MSH6wt",
                                   "K562"))

# categorize time variable into ordered factor (time_cat)
pData(eset) %<>%
  mutate(time_cat = ifelse(time == 0,
                           yes = "time_0",
                           no = ifelse(time %in% c(3,5,6,7,9,10),
                                       yes = "time_3_5_6_7_9_10",
                                       no = ifelse(time %in% c(12,14,15),
                                                   yes = "time_12_14_15",
                                                   # if including treated samples, there will be other time points, which should be accounted for
                                                   no = "Warning: time point other than 0,3,5,6,7,9,10,12,14,15"))),
         # set time_cat as ordered factor, and levels order
         time_cat = factor(time_cat,
                           ordered = TRUE,
                           levels = c("time_0",
                                      "time_3_5_6_7_9_10",
                                      "time_12_14_15")))

## list of conditions (column names in 'sampleinfo'):
# "cell_line","TP53","HMCES","A3A_vector","treatment","time","replicate","treatment_recat","time_cat"
conditions = names(pData(eset))[! names(pData(eset)) %in% c("sample", "time")]


### add library info to eset's 'featureData' section

# make sure order of grna_library's rows is the same as in the eset, i.e. 'grna_order'
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
#load(file = paste0(dataDir, "eset.RData")) ; conditions = names(pData(eset))[! names(pData(eset)) %in% c("sample", "time")]




#### NB regression

# IMPORTANT: depending on the analysis, change formula in src/functions.R

### run NB regression (loop across cell lines and genes) --> for full A3A dataset, it takes 8? hours in fsupeksvr

cell_line_names = pData(eset) %>%
  select(cell_line) %>%
  distinct() %>%
  filter(! str_detect(cell_line, "H358")) %>% # H358 does not have control samples
  pull(cell_line)

gene_names = fData(eset) %>%
  select(gene) %>%
  filter(gene != "Non-Targeting") %>%
  distinct() %>%
  pull(gene)

type_treatment = "control" # for the moment, only controls will be used

### load NB regression custom functions
source("src/functions.R")

NBres = list()

# turning off warning messages globally for the loop
options(warn=-1)

print(paste0("Running NB regressions on " , length(gene_names), " genes for ", length(gene_names), " cell lines:"))
             
for (cLine in cell_line_names){

  print(paste0("Analyzing cell line ", cLine))
  
  NBres_1_cell_line = mclapply(gene_names, function(g){
    
    # build dataframe (for a given gene) for regression
    df = buildDF(g, eset, conditions, raw_counts_offset, cLine, type_treatment)
    
    # run NB regression on that gene
    y = fitModel(g, eset, conditions, df)
    
  }, mc.cores = 5)
  
  # name elements by their analyzed gene's name
  names(NBres_1_cell_line) = gene_names
  
  # append to full NBres
  NBres[[cLine]] = NBres_1_cell_line
}

print("NB regressions finished!")

# warnings back
options(warn=0)

## save processed data
save(NBres, file = paste0(dataDir, "NBres.RData"))
#load(file = paste0(dataDir, "NBres.RData"))



# #### parse NB results
# 
# ## merge all genes results (i.e. the elements from NBres list) into megatable 'allres'
# 
# # apply broom::tidy() to each gene from NBres
# NBres_tidy = NBres %>%
#   map(., ~tidy(.x$model))
# 
# # now merge NB results for all genes into a megatable (takes time)
# allres = eat(NBres_tidy[[1]],
#              NBres_tidy[-1],
#              .by = "term") %>%
#   # the first dataset's variables remained unchanged, so append dataset name to them
#   rename_with(.fn = ~paste0(names(NBres_tidy)[1],
#                             "_estimate"),
#               .cols = all_of("estimate")) %>%
#   rename_with(.fn = ~paste0(names(NBres_tidy)[1],
#                             "_std.error"),
#               .cols = all_of("std.error")) %>%
#   rename_with(.fn = ~paste0(names(NBres_tidy)[1],
#                             "_statistic"),
#               .cols = all_of("statistic")) %>%
#   rename_with(.fn = ~paste0(names(NBres_tidy)[1],
#                             "_p.value"),
#               .cols = all_of("p.value")) %>%
#   # long format, by term
#   gather(key = "estimate.SE.stat.pval",
#          value = "value",
#          -term) %>%
#   # split 'estimate.SE.stat.pval' into 'estimate.SE.stat.pval' and 'gene' columns
#   separate(estimate.SE.stat.pval,
#            c("gene",
#              "estimate.SE.stat.pval"),
#            sep = "\\_") %>%
#   # one column estimate, one column SE, another statistic, another p.value
#   pivot_wider(names_from = c(estimate.SE.stat.pval),
#               values_from = c(value)) %>%
#   # remove NA rows
#   na.omit() %>%
#   # remove intercepts
#   filter(term != "(Intercept)") %>%
#   # add column with the condition tested (e.g. treatment_recatDOX_IC25-vs-control.. is "treatment") OR if it is interaction (contains ":")
#   mutate(condition = ifelse(str_detect(term, "\\:"),
#                             yes = "interaction",
#                             no = gsub("_.*",
#                                       "",
#                                       term))) %>%
#   # correct p-values (FDR)
#   mutate(fdr = p.adjust(p.value,
#                         method = "BH")) %>%
#   # order and sort
#   select(gene, condition, term, estimate, std.error, statistic, p.value, fdr) %>%
#   arrange(fdr, p.value)
# 
# ## write table to file
# write_tsv(allres,
#           paste0(tabDir, "results.tsv"))
# 
# 
# 
# #### detect top gene hits
# 
# # to detect genes with just negative slope for treated samples
# hits = data.frame(matrix(nrow=0,
#                          ncol=length(colnames(allres))))
# names(hits) = colnames(allres)
# hits = as_tibble(hits)
# for (g in unique(allres$gene)){
#   treatment_L = allres %>% filter(gene == g & term == 'treatment_recat.L')
#   time_mid_0 = allres %>% filter(gene == g & term == 'time_catmiddle-vs-zero')
#   time_late_mid = allres %>% filter(gene == g & term == 'time_catlate-vs-middle')
#   interaction_mid_0 = allres %>% filter(gene == g & term == 'time_catmiddle-vs-zero:treatment_recat.L')
#   interaction_late_mid = allres %>% filter(gene == g & term == 'time_catlate-vs-middle:treatment_recat.L')
#   ifelse(interaction_mid_0$fdr < 0.05 ||
#            interaction_late_mid$fdr < 0.05,
#          yes = (hits = rbind(hits, interaction_mid_0,
#                              interaction_late_mid,
#                              treatment_L,
#                              time_late_mid,
#                              time_mid_0)),
#          no = (hits = hits))
# }
# 
# 
# ## plot hits
# 
# for(gene in unique(hits$gene)){
#   data_df = NBres[[gene]]$df
#   data_model = NBres[[gene]]$model
#   new_model_names = c('counts', 'time_cat', 'treatment_recat', 'ln_sum_nontargeting')
#   names(data_model$model) = new_model_names
#   p_model = plot_model(data_model, 
#                        type = "pred", 
#                        terms = c('time_cat', 'treatment_recat'),
#                        value.offset = 'ln_sum_nontargeting') +
#     geom_line() +
#     ggtitle(gene) +
#     theme_minimal()
#   show(p_model)
#   ggsave(path = figDir, filename = paste0(gene, "_counts_vs_time.jpg"),
#           plot = p_model, device = "jpg", width = 10, height = 5.6, dpi = 300)  
# }

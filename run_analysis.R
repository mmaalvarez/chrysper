#### load libraries
library(vroom)
library(magrittr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)
library(readr)
library(DESeq2)
library(MASS)
library(ggplot2)
library(broom)
library(safejoin)
# library(effects)          # plot interaction terms
library(sjmisc)           # for plot_model function
library(sjPlot)           # for plot_model function

# resolve function conflicts
library(conflicted)
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")



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


# 350 non-essential genes (MAGeCK-defined) for offset
nonessentials = vroom("/g/strcombio/fsupek_cancer3/malvarez/genes_info/essentiality/core_sets/nonessential/mageck/mageck_nonessential_genes_count_normalization.txt", col_names = F, delim = " ")$X1



################################################################################################

### 1) Brunello datasets


### load datasets (tables of raw read counts and sample info), several options:

# raw data's root path
rawPath = "/g/strcombio/fsupek_data/CRISPR/"

### load sgRNAs library (Brunello), common for the 3 datasets
brunello_library = vroom(paste0(rawPath, "4_resources/crispr_libraries/brunello/original/broadgpp-brunello-library-contents_gRNAs.tsv")) %>%
  # update control sgRNAs name (non-targeting or non-essential, for offset later)
  mutate(gene = ifelse(str_detect(gene, "Non-Targeting") |
                         gene %in% nonessentials,
                       "control",
                       gene)) %>%
  # combine sgRNA and gene into an id
  unite(col = "sgRNA_id", sgRNA, gene, sep = "__", remove = F) %>%
  # keep interesting columns
  select(sgRNA_id, sgRNA, gene, "sgRNA Target Sequence", "Target Context Sequence", "PAM Sequence", "Strand", "Genomic Sequence", "Position of Base After Cut (1-based)") %>%
  # sort by chromosome and cut position
  arrange(`Genomic Sequence`, `Position of Base After Cut (1-based)`)

### raw data from 3 different experiments:

############################################################
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


####################################################
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


###############################################
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
  merge(select(brunello_library, sgRNA, `sgRNA Target Sequence`),
        by = "sgRNA Target Sequence") %>%
  select(-`sgRNA Target Sequence`) %>%
  relocate(sgRNA, .before = "gene")

# load sample info 
sampleinfo_oxygen = vroom(paste0(oxygen_Dir, "sampleinfo.tsv"))

# more info --> readLines(paste0(oxygen_Dir, "README"))


#### combine datasets
raw_counts_3_brunello_datasets = merge(raw_counts_A3A, raw_counts_ATRi_TMZ) %>%
  merge(raw_counts_oxygen) %>%
  # mark genes or non-targeting sgRNAs used later for offset
  mutate(gene = ifelse(gene == "Non-Targeting" |
                         gene %in% nonessentials,
                       "control",
                       gene))

sampleinfo_3_brunello_datasets = merge(sampleinfo_A3A, sampleinfo_ATRi_TMZ, all = T) %>%
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

sampleinfo_3_brunello_datasets %>%
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

## use the ln(sum(Non-Targeting sgRNA and 350 non-essential genes counts)) as an offset in the NB regression
# store this offset in a table for later adding it to the data frames created in buildDF()
raw_counts_offset_brunello = raw_counts_3_brunello_datasets %>%
  # long format
  gather(key = "sample",
         value = "counts",
         -c(sgRNA, gene)) %>%
  group_by(sample) %>%
  # create offset column
  mutate(ln_sum_control = log(sum(ifelse(gene == "control",
                                              yes = counts,
                                              no = 0),
                                       na.rm = TRUE))) %>%
  select(sample, ln_sum_control) %>%
  distinct()


#### store data as Biobase's 'ExpressionSets' (esets): a Class to Contain and Describe High-Throughput Expression Level Assays

# prepare raw counts table for conversion into eset
raw_counts_for_eset_brunello = raw_counts_3_brunello_datasets %>%
  # combine sgRNA and gene into an id
  unite(col = "sgRNA_id", sgRNA, gene, sep = "__") %>%
  # sgRNA_id to row names
  column_to_rownames("sgRNA_id") %>%
  # ExpressionSet requires a matrix object as input
  as.matrix()

eset_brunello = new("ExpressionSet", exprs = raw_counts_for_eset_brunello) # in case of wanting to use the vst normalized data, use 'exprs = vst'

# which is the order of sample names in eset
sample_order = colnames(eset_brunello)

# which is the order of sgRNA ids in eset
grna_order = rownames(eset_brunello) # takes the first column


### add sample info to eset's 'phenoData' section

# make sure order of sampleinfo's rows is the same as in the eset, i.e. 'sample_order'
# only samples that are BOTH in the counts AND sampleinfo tables are included (in theory, all of them)
sampleinfo_sorted = left_join(data.frame(sample = sample_order),
                              sampleinfo_3_brunello_datasets,
                              by = "sample")

pData(eset_brunello) = sampleinfo_sorted # columns: sample, id, cell_line, TP53, HMCES, A3A_vector, treatment, time, replicate


## set the phenoData's conditions that are categorical variables as factors: either unordered (use 'dummy' contrasts by default), or unordered (treatment and time, set contrasts below)
eset_brunello$cell_line = factor(eset_brunello$cell_line,
                        ordered = FALSE,
                        levels = c("A549_TP53KO_HMCESwt",
                                   "A549_TP53wt_HMCESwt",
                                   "LXF289_TP53mut_HMCESKO", # this cell line does not have control samples (woDOX + A3A, but this might be not a true control)
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
pData(eset_brunello) %<>%
  mutate(time_cat = ifelse(time == 0,
                           yes = "time_0",
                           no = ifelse(time %in% c(3:10),
                                       yes = "time_mid",
                                       no = ifelse(time %in% c(12:21),
                                                   yes = "time_late",
                                                   # if including treated samples, there will be other time points, which should be accounted for
                                                   no = "Warning: time point other than 0, 3:10, 12:21"))),
         # set time_cat as ordered factor, and levels order
         time_cat = factor(time_cat,
                           ordered = TRUE,
                           levels = c("time_0",
                                      "time_mid",
                                      "time_late")))

## list of conditions (column names in 'sampleinfo'):
conditions_brunello = names(pData(eset_brunello))[! names(pData(eset_brunello)) %in% c("sample", "time")]


### add library info to eset's 'featureData' section

# make sure order of grna_library's rows is the same as in the eset, i.e. 'grna_order'
# only sgRNAs that are BOTH in the counts AND library tables are included (in theory, all of them)
brunello_library_sorted = left_join(data.frame(sgRNA_id = grna_order),
                                brunello_library,
                                by = "sgRNA_id")

fData(eset_brunello) = brunello_library_sorted


## phenoData's 'sampleNames' section are just numbers -- replace with the actual sample names (SHOULD BE IN THE SAME ORDER)
sampleNames(eset_brunello) = eset_brunello$sample


## save processed data
save(eset_brunello,
     file = paste0(dataDir, "esets/eset_brunello.RData"))



#######################################################################################################################

#### 2) TKO v1 datasets

### load sgRNAs library (TKOv1), common for the 2 datasets
control_grnas = vroom("/g/strcombio/fsupek_data/CRISPR/4_resources/crispr_libraries/TKO/v1/control_sgRNAs_TKO_v1.txt", col_names = F, delim = " ") %>%
  separate(X1, into = c("prefix", "seq")) %>%
  pull(seq)
tkov1_library = vroom("/g/strcombio/fsupek_data/CRISPR/4_resources/crispr_libraries/TKO/v1/full/base/GSE128210_TKOv1-lib1-Human-Library.txt") %>%
  rename("gene" = "GENE",
         "sgRNA Target Sequence" = "SEQUENCE") %>%
  unite(col = "sgRNA", gene, `sgRNA Target Sequence`, remove = F) %>% 
  # update control sgRNAs name
  mutate(gene = ifelse(`sgRNA Target Sequence` %in% control_grnas |
                         gene %in% nonessentials,
                       "control",
                       gene)) %>%
  # combine sgRNA and gene into an id
  unite(col = "sgRNA_id", sgRNA, gene, sep = "__", remove = F) %>%
  # keep interesting columns
  select(sgRNA_id, sgRNA, gene, "sgRNA Target Sequence", EXON) %>%
  # sort
  arrange(gene, sgRNA)

## for final FDR table, which genes overlap between libraries
genes_overlap_libraries = merge(select(brunello_library, gene) %>% distinct(),
                                select(tkov1_library, gene) %>% distinct()) %>%
  pull(gene)


# raw data's root path
tkov1_rawPath = "/g/strcombio/fsupek_cancer3/malvarez/public_CRISPR_data/1_raw_data/"


#############################################
## iv) Moffat lab
# cell lines --> RPE1
# treatments --> control
# altered genotypes --> no (p53wt)

moffat_Dir = paste0(tkov1_rawPath, "isogenic_pairs_TP53/RPE1_moffat_lab/")

# load raw counts table
raw_counts_moffat = vroom(paste0(moffat_Dir, "GSE128210_Hart_et_al_Cell2015_RPE1-TKOv1-readcounts.txt")) %>%
  rename("sgRNA" = "GENE_CLONE",
         "gene" = "GENE") %>%
  separate(sgRNA, into = c("gene2", "seq"), sep = "_", remove = F) %>%
  # update control sgRNAs name
  mutate(gene = ifelse(seq %in% control_grnas |
                         gene %in% nonessentials,
                       "control",
                       gene)) %>%
  select(-c(seq, gene2))

# load sample info 
sampleinfo_moffat = vroom(paste0(moffat_Dir, "sample_info_RPE1_p53wt_tkov1_moffat.tsv"))

########################################
## v) Durocher lab
# cell lines --> RPE1, HeLa, SUM149PT
# treatments --> control
# altered genotypes --> no

durocher_Dir = paste0(tkov1_rawPath, "olaparib/DD0001_HeLa_RPE1_SUM149PT/")

# load raw counts table
raw_counts_durocher_RPE1 = vroom(paste0(durocher_Dir, "DD001_RPE1_CTRL.readcounts.txt")) %>%
  separate(GENE_CLONE, into = c("gene", "seq"), sep = "_", remove = F) %>%
  # update control sgRNAs name
  mutate(gene = ifelse(seq %in% control_grnas |
                         gene %in% nonessentials,
                       "control",
                       gene)) %>%
  select(-seq) %>%
  rename("sgRNA" = "GENE_CLONE")

raw_counts_durocher_HeLa = vroom(paste0(durocher_Dir, "DD001_HeLa_CTRL.withrepBT15.readcounts.txt")) %>%
  separate(GENE_CLONE, into = c("gene", "seq"), sep = "_", remove = F) %>%
  # update control sgRNAs name
  mutate(gene = ifelse(seq %in% control_grnas |
                         gene %in% nonessentials,
                       "control",
                       gene)) %>%
  select(-seq) %>%
  rename("sgRNA" = "GENE_CLONE")

raw_counts_durocher_SUM149PT = vroom(paste0(durocher_Dir, "DD001_SUM149PT_CTRL.readcounts.txt")) %>%
  separate(GENE_CLONE, into = c("gene", "seq"), sep = "_", remove = F) %>%
  # update control sgRNAs name
  mutate(gene = ifelse(seq %in% control_grnas |
                         gene %in% nonessentials,
                       "control",
                       gene)) %>%
  select(-seq) %>%
  rename("sgRNA" = "GENE_CLONE")

# load sample info 
sampleinfo_durocher = vroom(paste0(durocher_Dir, "sample_info.tsv"))


#### combine datasets
raw_counts_2_tkov1_datasets = merge(raw_counts_moffat, raw_counts_durocher_RPE1) %>%
  merge(raw_counts_durocher_HeLa) %>%
  merge(raw_counts_durocher_SUM149PT)

sampleinfo_2_tkov1_datasets = merge(sampleinfo_moffat, sampleinfo_durocher, all = T) %>%
  # control or treated
  mutate(type = ifelse(treatment == "no",
                       "control",
                       "treated")) %>%
  select(sample, cell_line, time, type, replicate) %>%
  arrange(cell_line, time, type, replicate)

sampleinfo_2_tkov1_datasets %>%
  filter(type == "control") %>%
  select(cell_line, time) %>%
  table()
#     time 0 9 12 15 18 21
#---------------------------
# HeLa     1 3  0  3  0  3
# RPE1     2 4  4  4  4  2
# SUM149PT 1 3  3  3  3  3


#### normalize count matrix (to allow comparisons between samples)

## use the ln(sum(control sgRNA and 350 non-essential genes counts)) as an offset in the NB regression
# store this offset in a table for later adding it to the data frames created in buildDF()
raw_counts_offset_tkov1 = raw_counts_2_tkov1_datasets %>%
  # long format
  gather(key = "sample",
         value = "counts",
         -c(sgRNA, gene)) %>%
  group_by(sample) %>%
  # create offset column
  mutate(ln_sum_control = log(sum(ifelse(gene == "control",
                                              yes = counts,
                                              no = 0),
                                       na.rm = TRUE))) %>%
  select(sample, ln_sum_control) %>%
  distinct()



#### store data as Biobase's 'ExpressionSets' (esets)

# prepare raw counts table for conversion into eset
raw_counts_for_eset_tkov1 = raw_counts_2_tkov1_datasets %>%
  # combine sgRNA and gene into an id
  unite(col = "sgRNA_id", sgRNA, gene, sep = "__") %>%
  # sgRNA_id to row names
  column_to_rownames("sgRNA_id") %>%
  # ExpressionSet requires a matrix object as input
  as.matrix()

eset_tkov1 = new("ExpressionSet", exprs = raw_counts_for_eset_tkov1) # in case of wanting to use the vst normalized data, use 'exprs = vst'

# which is the order of sample names in eset
sample_order = colnames(eset_tkov1)

# which is the order of sgRNA ids in eset
grna_order = rownames(eset_tkov1) # takes the first column


### add sample info to eset's 'phenoData' section

# make sure order of sampleinfo's rows is the same as in the eset, i.e. 'sample_order'
# only samples that are BOTH in the counts AND sampleinfo tables are included (in theory, all of them)
sampleinfo_sorted = left_join(data.frame(sample = sample_order),
                              sampleinfo_2_tkov1_datasets,
                              by = "sample")

pData(eset_tkov1) = sampleinfo_sorted # columns: sample, id, cell_line, TP53, HMCES, A3A_vector, treatment, time, replicate


## set the phenoData's conditions that are categorical variables as factors: either unordered (use 'dummy' contrasts by default), or unordered (treatment and time, set contrasts below)
eset_tkov1$cell_line = factor(eset_tkov1$cell_line,
                                 ordered = FALSE,
                                 levels = c("RPE1",
                                            "HeLa",
                                            "SUM149PT"))

# categorize time variable into ordered factor (time_cat)
pData(eset_tkov1) %<>%
  mutate(time_cat = ifelse(time == 0,
                           yes = "time_0",
                           no = ifelse(time %in% c(3:10),
                                       yes = "time_mid",
                                       no = ifelse(time %in% c(12:21),
                                                   yes = "time_late",
                                                   # if including treated samples, there will be other time points, which should be accounted for
                                                   no = "Warning: time point other than 0, 3:10, 12:21"))),
         # set time_cat as ordered factor, and levels order
         time_cat = factor(time_cat,
                           ordered = TRUE,
                           levels = c("time_0",
                                      "time_mid",
                                      "time_late")))

## list of conditions (column names in 'sampleinfo'):
conditions_tkov1 = names(pData(eset_tkov1))[! names(pData(eset_tkov1)) %in% c("sample", "time")]


### add library info to eset's 'featureData' section

# make sure order of grna_library's rows is the same as in the eset, i.e. 'grna_order'
# only sgRNAs that are BOTH in the counts AND library tables are included (in theory, all of them)
tkov1_library_sorted = left_join(data.frame(sgRNA_id = grna_order),
                                    tkov1_library,
                                    by = "sgRNA_id")

fData(eset_tkov1) = tkov1_library_sorted

## phenoData's 'sampleNames' section are just numbers -- replace with the actual sample names (SHOULD BE IN THE SAME ORDER)
sampleNames(eset_tkov1) = eset_tkov1$sample

## save processed data
save(eset_tkov1,
     file = paste0(dataDir, "esets/eset_tkov1.RData"))




#######################################################################################################################


#### Prepare NB regressions

# load eset data
load(file = paste0(dataDir, "esets/eset_brunello.RData")) ; conditions_brunello = names(pData(eset_brunello))[! names(pData(eset_brunello)) %in% c("sample", "time")]
load(file = paste0(dataDir, "esets/eset_tkov1.RData")) ; conditions_tkov1 = names(pData(eset_tkov1))[! names(pData(eset_tkov1)) %in% c("sample", "time")]

# cell line and gene names per library
cell_line_names_brunello = pData(eset_brunello) %>%
  select(cell_line) %>%
  distinct() %>%
  filter(! str_detect(cell_line, "H358|LXF289_TP53mut_HMCESKO")) %>% # H358 and LXF289_TP53mut_HMCESKO do not have control samples
  pull(cell_line)
cell_line_names_tkov1 = pData(eset_tkov1) %>%
  select(cell_line) %>%
  distinct() %>%
  pull(cell_line)
gene_names_brunello = fData(eset_brunello) %>%
  select(gene) %>%
  filter(gene != "control") %>%
  distinct() %>%
  pull(gene)
gene_names_tkov1 = fData(eset_tkov1) %>%
  select(gene) %>%
  filter(gene != "control") %>%
  distinct() %>%
  pull(gene)


#### Start NB regressions

### load NB regression custom functions
source("src/functions.R")

### start with Brunello data
print("Starting with Brunello datasets")

## first do regressions for each gene, with all (control) cell lines pooled
print(paste0("Running NB regressions on " , length(gene_names_brunello), " genes for all ", length(cell_line_names_brunello), " cell lines pooled"))

# NB loop across genes
NBres = list()
NBres[["Brunello_pooled"]] = lapply(gene_names_brunello, function(g){ #mc
    
  # build dataframe (for a given gene) for regression
  df = buildDF(g, eset_brunello, conditions_brunello, raw_counts_offset_brunello, cLine="pooled", type_treatment="control")
  
  # run NB regression on that gene
  y = fitModel(df)
    
}) #, mc.preschedule=FALSE, mc.cores = 8) # with preschedule=T, >1 cores gives empty values (errors?) for some genes...
  
# name elements by their analyzed gene's name
names(NBres[["Brunello_pooled"]]) = gene_names_brunello


## now do regressions for each gene, for each cell line separately
print(paste0("Now running NB regressions on " , length(gene_names_brunello), " genes for ", length(cell_line_names_brunello), " cell lines separately"))

# NB loop across cell lines and genes
for (cLine in cell_line_names_brunello){

  print(paste0("Analyzing cell line ", cLine))
  
  NBres_1_cell_line = lapply(gene_names_brunello, function(g){ #mc
    
    # build dataframe (for a given gene) for regression
    df = buildDF(g, eset_brunello, conditions_brunello, raw_counts_offset_brunello, cLine, type_treatment="control")
    
    # run NB regression on that gene
    y = fitModel(df)
    
  }) #, mc.preschedule=FALSE, mc.cores = 8)
  
  # name elements by their analyzed gene's name
  names(NBres_1_cell_line) = gene_names_brunello
  
  # append to full NBres
  NBres[[cLine]] = NBres_1_cell_line
  
  # free memory
  rm(NBres_1_cell_line); gc()
}


### now TKO v1
print("Continue with TKO v1 datasets")

## first do regressions for each gene, with all (control) cell lines pooled
print(paste0("Running NB regressions on " , length(gene_names_tkov1), " genes for all ", length(cell_line_names_tkov1), " cell lines pooled"))

# NB loop across genes
NBres[["TKOv1_pooled"]] = lapply(gene_names_tkov1, function(g){ #mc
  
  # build dataframe (for a given gene) for regression
  df = buildDF(g, eset_tkov1, conditions_tkov1, raw_counts_offset_tkov1, cLine="pooled", type_treatment="control")
  
  # run NB regression on that gene
  y = fitModel(df)
  
}) #, mc.preschedule=FALSE, mc.cores = 8)

# name elements by their analyzed gene's name
names(NBres[["TKOv1_pooled"]]) = gene_names_tkov1


## now do regressions for each gene, for each cell line separately
print(paste0("Running NB regressions on " , length(gene_names_tkov1), " genes for ", length(cell_line_names_tkov1), " cell lines separately"))

# NB loop across cell lines and genes
for (cLine in cell_line_names_tkov1){
  
  print(paste0("Analyzing cell line ", cLine))
  
  NBres_1_cell_line = lapply(gene_names_tkov1, function(g){ #mc
    
    # build dataframe (for a given gene) for regression
    df = buildDF(g, eset_tkov1, conditions_tkov1, raw_counts_offset_tkov1, cLine, type_treatment="control")
    
    # run NB regression on that gene
    y = fitModel(df)
    
  }) #, mc.preschedule=FALSE, mc.cores = 8)
  
  # name elements by their analyzed gene's name
  names(NBres_1_cell_line) = gene_names_tkov1
  
  # append to full NBres
  NBres[[cLine]] = NBres_1_cell_line
  
  # free memory
  rm(NBres_1_cell_line); gc()
}

print("NB regressions finished!")

# ## save processed data
# save(NBres, file = paste0(dataDir, "NBres/NBres.RData"),
#      compress = F) # the parameter compress=T means that the resulting file will use less space on your disk. However, if it is a really huge dataset, it could take longer to load it later because R first has to extract the file again. So, if you want to save space, then leave it as it is. If you want to save time, add a parameter compress = F

# NBres size
# 13G in disk (dataDir/NBres/)
# 25G in memory (workspace)

# free memory up
rm(list = ls(pattern="grna|library|sampleinfo|eset|raw|gene|line|table"))
gc()



#### parse NB results

# load processed data
load(file = paste0(dataDir, "NBres/NBres.RData"))

## loop through NB list
for (cLine in names(NBres)){
  
  NBres_tidy = NBres[[cLine]] %>%
    # apply broom::tidy() to each gene from NBres[[cell_line]] -- it's like summary..
    map(., ~tidy(.x)) %>%
    # keeping the estimate and pval of time_cat.Q only
    map(., ~filter(.x, term == "time_cat.Q")) %>%
    map(., ~select(.x, -c(term, std.error, statistic))) %>%
    # add gene name
    map2_df(., names(.), ~mutate(.x, gene = .y)) %>%
    relocate(gene, .before="estimate") %>%
    # add cell line name
    mutate(cell_line = cLine) %>%
    relocate(cell_line, .after="gene")
  
  # write to intermediate tables
  write_tsv(NBres_tidy, paste0(dataDir, "estimates_pvalues/estimates_pvalues_", cLine, ".tsv" ))

  # free memory
  rm(NBres_tidy); gc()
}

# free memory
rm(NBres); gc()

## recover tables for each cell line and rbind them
list_tables = lapply(Sys.glob(paste0(dataDir, "estimates_pvalues/estimates_pvalues_*")),
                                     vroom)
table_fdr = eat(list_tables[[1]],
                list_tables[-1],
                .mode = "full") %>%
  # correct p.values (FDR)
  mutate(fdr = qvalue::qvalue(p.value)$qvalues, #p.adjust(p.value, method = "BH",
         # mark which genes do not overlap between Brunello and TKOv1)
         library_overlap = ifelse(gene %in% genes_overlap_libraries,
                                  "Library overlap",
                                  "No library overlap"),
         # count trend shape, for faceting
         shape = ifelse(estimate>0,
                               "convex",
                               "concave"),
         shape = factor(shape, ordered = T, levels = c("convex", "concave"))) %>%
  # sort based on gene
  arrange(gene)

# free memory
rm(list_tables); gc()

# write table to file
write_tsv(table_fdr, paste0(tabDir, "table_fdr.tsv"))



## gene hit ascertainment (FDR<..., overlapping across many control cell lines)

# load table_fdr
table_fdr = vroom(paste0(tabDir, "table_fdr.tsv"))

# threshold is based on the mean(-log10(fdr)) across cell lines, per gene AND SHAPE
minlog10fdr_thr = 0.60206 # fdr = 25% (-log10(0.25)=0.60206)

# point size = mean(-log10(fdr)) across cell lines, per gene AND SHAPE
minlog10fdr_mean = table_fdr %>%
  select(gene, shape, fdr) %>%
  group_by(gene, shape) %>%
  summarize(fdr_mean = mean(fdr),
            .groups = "keep") %>%
  mutate(minlog10fdr_mean = -log10(fdr_mean)) %>%
  select(-fdr_mean)
table_fdr %<>% merge(minlog10fdr_mean)

# plot
manhattan = ggplot(data = table_fdr,
       aes(x = gene,
           # -log scale
           y = -log10(fdr))) +
  scale_y_continuous(breaks = c(seq(0, 300, 25)),
                     expand = c(0, 0)) +
  geom_point(alpha = 0.8,
             aes(col = cell_line,
                 shape = library_overlap,
                 size = minlog10fdr_mean)) +
  scale_shape_manual(values = c(16, 18)) +
  geom_text_repel(data = table_fdr %>% #filter(label_this == "label_this"),
                    filter(minlog10fdr_mean >= minlog10fdr_thr &
                             -log10(fdr) >= minlog10fdr_thr),
                  aes(label = gene),
                  nudge_y = 150,
                  force = 2,
                  size = 2,
                  min.segment.length = 0.1,
                  segment.size = 0.1,
                  show.legend = F) +
  facet_wrap(facets = ~shape,
             scales = "free") +
  xlab("gene") +
  ylab("-log10(q-value)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5,'cm'), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=5)) #change legend text font size
ggsave(path = figDir, plot = manhattan,
       filename = "manhattan.jpg",
       device = "jpg", width = 10, height = 5.6, dpi = 600)
ggsave(path = figDir, plot = manhattan,
       filename = "manhattan.pdf",
       device = "pdf", width = 10, height = 5.6, dpi = 600)

manhattan_clean = ggplot(data = table_fdr %>%
                     # plot only genes~cLine with FDR<thr
                     filter(minlog10fdr_mean >= minlog10fdr_thr),
                   aes(x = gene,
                       # -log scale
                       y = -log10(fdr))) +
  scale_y_continuous(breaks = c(seq(0, 300, 25)),
                     expand = c(0, 0)) +
  geom_point(alpha = 0.8,
             aes(col = cell_line,
                 shape = library_overlap,
                 size = minlog10fdr_mean)) +
  scale_shape_manual(values = c(16, 18)) +
  geom_text_repel(data = table_fdr %>% #filter(label_this == "label_this"),
                    filter(minlog10fdr_mean >= minlog10fdr_thr &
                             -log10(fdr) >= minlog10fdr_thr),
                  aes(label = gene),
                  nudge_y = 50,
                  force = 2,
                  size = 2,
                  min.segment.length = 0.1,
                  segment.size = 0.1,
                  show.legend = F) +
  facet_wrap(facets = ~shape,
             scales = "free") +
  xlab("gene") +
  ylab("-log10(q-value)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5,'cm'), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=5)) #change legend text font size
ggsave(path = figDir, plot = manhattan_clean,
       filename = "manhattan_clean.jpg",
       device = "jpg", width = 10, height = 5.6, dpi = 600)
ggsave(path = figDir, plot = manhattan_clean,
       filename = "manhattan_clean.pdf",
       device = "pdf", width = 10, height = 5.6, dpi = 600)


## draw hit curves with plot_model

load(file = paste0(dataDir, "NBres/NBres.RData"))

# the hits are the gene~cLine~shape shown in manhattan_clean
table_fdr$shape = factor(table_fdr$shape, ordered = T, levels = c("convex", "concave"))
hits = table_fdr %>%
  filter(minlog10fdr_mean >= minlog10fdr_thr &
           -log10(fdr) >= minlog10fdr_thr) %>%
  select(gene, cell_line, shape, library_overlap) %>%
  arrange(shape, gene, cell_line)

list_plots = list()
for(i in seq(1, length(hits$gene))){
  data_model = NBres[[hits$cell_line[i]]][[hits$gene[i]]]
  # new model col names
  names(data_model$model) = c('counts', 'time_cat', 'ln_sum_control')
  p_model = plot_model(data_model,
                       type = "pred",
                       terms = c('time_cat'),
                       value.offset = 'ln_sum_control') +
    xlab("") +
    geom_line() +
    ggtitle(paste0(hits$gene[i], " (", hits$library_overlap[i], ") - ", hits$cell_line[i], " - ", hits$shape[i] , " shape")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  list_plots[[paste0(hits$gene[i], "_", hits$cell_line[i], "_", hits$shape[i])]] = p_model
}
#show(list_plots[["WDR52_RPE1_concave"]])

cowplot(list_plots)


rm(NBres); gc()



## GO set enrichment of gene hits can provide more insights


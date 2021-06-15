### create testing dataset

## using some A3A 1st experiment samples:
# sample            id    cell_line TP53  HMCES A3A_vector treatment  time  replicate
# 1 1A_p53KO_A3A_t0 1A    A549      KO    wt    yes        no            0         1
# 2 4A_p53KO_A3A_w… 4A    A549      KO    wt    yes        no            9         1
# 3 6A_p53KO_A3A_w… 6A    A549      KO    wt    yes        no           15         1
# 4 9A_p53KO_A3A_D… 9A    A549      KO    wt    yes        DOX_IC25      9         1
# 5 11A_p53KO_A3A_… 11A   A549      KO    wt    yes        DOX_IC25     15         1
# 6 14A_p53KO_A3A_… 14A   A549      KO    wt    yes        DOX_ATRi      9         1
# 7 16A_p53KO_A3A_… 16A   A549      KO    wt    yes        DOX_ATRi     15         1

## and 2 "random" genes + Non-Targeting sgRNAs:
genes = c("A1BG", "A1CF", "Non-Targeting")


## path
rawPath = "/g/strcombio/fsupek_data/CRISPR/"
A3A_Dir = paste0(rawPath, "2_processed_data/apobec/MAGeCK_MLE_analyses/raw_counts_combined/")


## load raw counts table and create test
A3A_raw_counts = vroom(paste0(A3A_Dir, "A549_LXF289_H358__raw_counts.tsv")) %>%
  filter(gene %in% genes) %>%
  select(sgRNA, gene, matches("[0-9]A_p53KO")) %>%
  select(sgRNA, gene, matches("A3A")) %>%
  select(sgRNA, gene, matches(c("t0", "t9", "t15")))

write_tsv(A3A_raw_counts,
          "test_raw_counts.tsv")


## load sample info and create test
A3A_sampleinfo = vroom(paste0(A3A_Dir, "sampleinfo.tsv")) %>%
  filter(sample %in% names(A3A_raw_counts))

write_tsv(A3A_sampleinfo,
          "test_sampleinfo.tsv")


### load test gene essentiality (mean DEMETER2 score <0 --> gene inhibition decreases fitness consistently across cell lines)
D2_score = vroom(paste0(rawPath, "4_resources/genes_info/essentiality/depmap/demeter/mean_D2.tsv")) %>%
  filter(gene %in% genes)


## load test sgRNA library
grna_library = vroom(paste0(rawPath, "4_resources/crispr_libraries/brunello/original/broadgpp-brunello-library-contents_gRNAs.tsv")) %>%
  rename("gene" = "Target Gene Symbol") %>%
  # update Non-Targeting sgRNAs name (remove trailing " Control" part)
  mutate(gene = gsub(" Control$", "", gene)) %>%
  # combine sgRNA and gene into an id
  unite(col = "sgRNA_id", sgRNA, gene, sep = "__", remove = F) %>%
  # keep interesting columns
  select(sgRNA_id, sgRNA, gene, "sgRNA Target Sequence", "Target Context Sequence", "PAM Sequence", "Strand", "Genomic Sequence", "Position of Base After Cut (1-based)") %>%
  # sort by chromosome and cut position
  arrange(`Genomic Sequence`, `Position of Base After Cut (1-based)`) %>%
  # extract 2 genes + non-targeting sgRNAs
  filter(gene %in% genes) %>%
  # add demeter2 score
  merge(D2_score, all = T)

write_tsv(grna_library,
          "test_grna_library.tsv")

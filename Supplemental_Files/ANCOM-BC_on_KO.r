library(ANCOMBC)
library(tidyverse)
#library(caret)
#library(DT)

library(phyloseq)
metadata <- read.delim("C:/Gresham/tiny_projects/Project_Osun/metadata/MappingFile_Osun_ANCOM.tsv")
counts_table <- read.delim("C:/Gresham/tiny_projects/Project_Osun/picrust/ko_picrust2_feature-table.tsv", row.names=1)

metadata<-metadata
metadata
str(metadata)
info_data<-counts_table[c(4, 5, 6, 7, 8)]
info_data
X <-phyloseq(
  otu_table(info_data, taxa_are_rows = T), 
  sample_data(metadata %>% as.data.frame() %>% column_to_rownames("sample.id")))
X  

set.seed(123)
output = ancombc2(data = X, assay_name = "counts", tax_level = NULL,
                  fix_formula = "Site", rand_formula = NULL,
                  p_adj_method = "fdr", pseudo = 0, pseudo_sens = TRUE,
                  prv_cut = 0.15, lib_cut = 1000, s0_perc = 0.05,
                  group = "Site", struc_zero = TRUE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = FALSE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, mdfdr_control = NULL, 
                  trend_control = NULL)

results <- output$res
write.table(results, file = "C:/Gresham/tiny_projects/Project_Osun/picrust/ANCOM_KO_Valley_v_Spring.txt", sep = "\t")

info_data<-counts_table[c(2, 3, 5, 6)]
#comp <- counts_table[c(1, 2, 4, 5)]
#info_data
X <-phyloseq(
  otu_table(info_data, taxa_are_rows = T), 
  sample_data(metadata %>% as.data.frame() %>% column_to_rownames("name")))
X  

set.seed(123)
output = ancombc2(data = X, assay_name = "counts", tax_level = NULL,
                  fix_formula = "exp", rand_formula = NULL,
                  p_adj_method = "fdr", pseudo = 0, pseudo_sens = TRUE,
                  prv_cut = 0.9, lib_cut = 1000, s0_perc = 0.05,
                  group = "exp", struc_zero = TRUE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = FALSE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, mdfdr_control = NULL, 
                  trend_control = NULL)

results <- output$res
write.table(results, file = "ANCOM_Obs_DGY1657_DGY1728_reps_23_v2.txt", sep = "\t")
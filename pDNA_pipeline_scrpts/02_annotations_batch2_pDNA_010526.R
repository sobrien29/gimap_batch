user_lib <- "~/R/library"
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)

# prepend it to R library paths
.libPaths(c(user_lib, .libPaths()))

# helper function to install missing packages
install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, lib = user_lib, repos = "https://cloud.r-project.org")
      library(pkg, character.only = TRUE)
    }
  }
}

# CRAN packages
cran_packages <- c("tidyverse", "tidylog", "kableExtra", "pheatmap", "RColorBrewer")
install_if_missing(cran_packages)

# Bioconductor packages
bioc_packages <- c("depmap", "biomaRt", "ExperimentHub")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = user_lib, repos = "https://cloud.r-project.org")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE, lib = user_lib)
    library(pkg, character.only = TRUE)
  }
}

message("All packages loaded successfully.")

library(tidyverse)
library(tidylog)
library("depmap") # for CN data
library(biomaRt)
library("ExperimentHub") 
library(kableExtra) # for formatting kables



args <- commandArgs(trailingOnly = TRUE)

base_dir  <- args[1]
cell_line <- args[2]



source(file.path(base_dir, "00-shared_functions_and_variables.R"))

in_dir <- file.path(base_dir, "config")
if (!dir.exists(in_dir)) dir.create(in_dir, recursive = TRUE)

out_dir <- file.path(base_dir, "results", "pgRNA_annotations")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


d.annot <- read_tsv(file.path(in_dir, "pgMI_annotations.txt"))
d.pan_essential <- read_csv(file.path(in_dir, "Achilles_common_essentials.csv"))
d.gene_name_to_id <- read_tsv(file.path(in_dir, "hgnc_to_ensembl.txt"))

## get just gene names and Ensembl IDs
d.genes <- d.annot %>%
  filter(target_type == "gene_gene") %>%
  dplyr::select(paralog_pair, paralog_pair_id) %>%
  separate_rows(paralog_pair, paralog_pair_id, sep = "_") %>%
  rename(gene_symbol = paralog_pair, ensembl_id = paralog_pair_id) %>%
  distinct(gene_symbol, .keep_all = TRUE)

## get vector of gene IDs
ensembl_ids <- d.genes %>%
  pull(ensembl_id)

ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=101) ## Alice had to update

## checking Ensembl settings
filters <- listFilters(ensembl)

## get paralog Entrez IDs

d.all_ids <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
                   filters = c("ensembl_gene_id"), 
                   values = list(ensembl_ids), 
                   mart = ensembl)

## which genes have multiple Entrez ID entries?
##   checking this by hand shows that all larger numbers are ncRNAs... 
##   so for now I will remove them using distinct()
d.all_ids %>%
  group_by(hgnc_symbol) %>%
  summarize(n = n()) %>%
  filter(n > 1)

## arrange by entrez_id within ensembl_id groups (will only affect genes w/ 1 Ensembl
##   but multiple Entrez IDs), then use distinct to keep just the first row
d.all_ids_unique <- d.all_ids %>%
  rename("entrez_id" = entrezgene_id, "ensembl_id" = ensembl_gene_id) %>%
  group_by(ensembl_id) %>%
  arrange(entrez_id, .by_group = TRUE) %>%
  distinct(ensembl_id, .keep_all = TRUE)
d.all_ids_unique
nrow(d.all_ids_unique)

## warning message
if(nrow(d.all_ids_unique) < nrow(d.genes)){
  warning("Warning: not all of your genes have matching Entrez IDs!")
}

## join Entrez and HGNC IDs with the original library gene list
d.pgMI_gene_ids_all <- d.genes %>%
  left_join(d.all_ids_unique, by = "ensembl_id") 

## get TPM and CN information (w/ option for user to upload their own info)
## add if statement for user TPM/CN vs. DepMap
## print label/version of dataset - save a tbl? 
## or use a specific version... hmm... 
d.depmap_tpm <- depmap::depmap_TPM()

d.depmap_metadata <- depmap::depmap_metadata()
d.depmap_metadata <- d.depmap_metadata %>% dplyr::select(-cell_line)

## store DepMap and EH IDs for later use
depmap_release_id <- depmap::depmap_release()
eh <- ExperimentHub()
eh_id <- names(query(eh, c("depmap", paste("TPM", depmap_release_id, sep = "_"))))

# ===============================
# DepMap & Annotation Processing
# ===============================

library(dplyr)
library(stringr)

# Ensure cell_line is a clean string
cell_line <- trimws(cell_line)
if(length(cell_line) != 1) stop("cell_line must be a single string")

# Find DepMap ID
my_depmap_id <- d.depmap_metadata %>%
  filter(str_detect(stripped_cell_line_name, fixed(cell_line, ignore_case = TRUE))) %>%
  pull(depmap_id)

# Flag whether DepMap data exists
depmap_available <- length(my_depmap_id) == 1

if(depmap_available){
  
  message("DepMap ID found for ", cell_line, ": ", my_depmap_id)
  
  # --- TPM ---
  d.depmap_tpm_my_cell_line <- d.depmap_tpm %>%
    filter(depmap_id == my_depmap_id) %>%
    dplyr::select(gene_name, entrez_id, rna_expression) %>%
    rename(depmap_gene_symbol = gene_name, log2_tpm = rna_expression)
  
  d.pgMI_gene_tpm <- d.pgMI_gene_ids_all %>%
    left_join(d.depmap_tpm_my_cell_line, by = "entrez_id")
  
  rm(d.depmap_tpm_my_cell_line)
  
  d.pgMI_gene_tpm_filtered <- d.pgMI_gene_tpm %>%
    group_by(gene_symbol) %>%
    mutate(duplicated = n() > 1) %>%
    mutate(remove = case_when(
      duplicated & gene_symbol != depmap_gene_symbol ~ TRUE,
      TRUE ~ FALSE
    )) %>%
    filter(remove == FALSE) %>%
    dplyr::select(-c(remove, duplicated)) %>%
    mutate(expressed_flag = case_when(
      log2_tpm < 1 ~ FALSE,
      log2_tpm >= 1 ~ TRUE,
      is.na(log2_tpm) ~ NA
    ))
  
  # TPM plot
  library_gene_tpm_with_cutoff <- ggplot(d.pgMI_gene_tpm_filtered, aes(x = log2_tpm)) +
    geom_histogram(binwidth = 0.5, color = "black", fill = "darkgray") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    plot_options + plot_theme + theme(aspect.ratio = wide_ar)
  save_plot(library_gene_tpm_with_cutoff)
  
  # --- CN ---
  d.depmap_cn_my_cell_line <- depmap::depmap_copyNumber() %>%
    filter(depmap_id == my_depmap_id) %>%
    dplyr::select(gene_name, entrez_id, log_copy_number) %>%
    rename(depmap_gene_symbol = gene_name, log2_cn = log_copy_number)
  
  d.pgMI_gene_cn <- d.pgMI_gene_tpm_filtered %>%
    left_join(d.depmap_cn_my_cell_line, by = "entrez_id")
  
  # CN plot
  library_gene_cn <- ggplot(d.pgMI_gene_cn, aes(x = log2_cn)) +
    geom_histogram(binwidth = 0.1, color = "black", fill = "darkgray") +
    plot_options + plot_theme + theme(aspect.ratio = wide_ar)
  save_plot(library_gene_cn)
  
  # Prepare df for annotation joins
  d.pgMI_gene_tpm_cn <- d.pgMI_gene_cn %>%
    dplyr::select(-c(depmap_gene_symbol.x, depmap_gene_symbol.y))
  
} else {
  
  if(length(my_depmap_id) == 0){
    message("No DepMap ID found for ", cell_line, ". Skipping TPM/CN processing.")
  } else {
    warning("Multiple DepMap IDs found for ", cell_line, ". Skipping TPM/CN processing.")
  }
  
  # Create placeholders so annotation joins still work
  d.pgMI_gene_tpm_filtered <- d.pgMI_gene_ids_all %>%
    mutate(log2_tpm = NA, expressed_flag = NA)
  d.pgMI_gene_tpm_cn <- d.pgMI_gene_tpm_filtered %>%
    mutate(log2_cn = NA)
}
## split up gene col into symbol and entrez ID
d.pan_essential <- d.pan_essential %>%
  separate(col = gene, into = c("gene_symbol", "entrez_id"), remove = FALSE)

pan_essential_entrez <- d.pan_essential %>%
  pull(entrez_id)

d.pgMI_gene_tpm_cn_ctrl <- d.pgMI_gene_tpm_cn %>%
  mutate(essential_flag = case_when(
    entrez_id %in% pan_essential_entrez ~ TRUE,
    TRUE ~ FALSE)) ## if entrez_id is NA or not in pan-essential Entrez ID list, FALSE 

## number of pan-essential genes in my library
d.pgMI_essential_summary <- d.pgMI_gene_tpm_cn_ctrl %>%
  group_by(essential_flag) %>%
  summarize(n = n())
print_kbl(d.pgMI_essential_summary)

d.annot <- d.annot %>%
  left_join(
    d.pgMI_gene_tpm_cn_ctrl %>%
      dplyr::select(ensembl_id, entrez_id:essential_flag) %>%
      rename_with(~ paste("gene1", .x, sep = "_")),
    by = c("gene1_ensembl_id")
  ) %>%
  left_join(
    d.pgMI_gene_tpm_cn_ctrl %>%
      dplyr::select(ensembl_id, entrez_id:essential_flag) %>%
      rename_with(~ paste("gene2", .x, sep = "_")),
    by = c("gene2_ensembl_id")
  )

d.annot <- d.annot %>%
  distinct(pgRNA_id, gene1_ensembl_id, gene2_ensembl_id, .keep_all = TRUE)

d.annot

d.annot <- d.annot %>%
  mutate(norm_ctrl_flag = case_when(
    target_type == "gene_gene" ~ "double_targeting",
    target_type == "gene_ctrl" & gene1_essential_flag == TRUE ~ "positive_control",
    target_type == "ctrl_gene" & gene2_essential_flag == TRUE ~ "positive_control",
    target_type == "gene_ctrl" & gene1_essential_flag != TRUE ~ "single_targeting", 
    target_type == "ctrl_gene" & gene2_essential_flag != TRUE ~ "single_targeting",
    target_type == "ctrl_ctrl" ~ "negative_control")) %>%
  mutate(norm_ctrl_flag = factor(norm_ctrl_flag, levels = c("negative_control",
                                                            "positive_control",
                                                            "single_targeting",
                                                            "double_targeting")))

## add a flag 
d.annot_norm_ctrl_summary <- d.annot %>%
  group_by(norm_ctrl_flag) %>%
  summarize(n = n())
print_kbl(d.annot_norm_ctrl_summary)

d.annot_neg_ctrl_only <- d.annot %>%
  filter(target_type == "ctrl_ctrl") %>%
  dplyr::select(pgRNA_id)

## making fake "genes" for ctrl_ctrl pgRNAs so I can normalize like the others
## Q: is this necessary if I'm not using MAGeCK? 
set.seed(123)
N <- nrow(d.annot_neg_ctrl_only)
fake_gene_nums <- rep(seq(from = 1, to = 50, by = 1), 10)

# randomize the vector
fake_gene_nums <- sample(rep(seq(1, 50), length.out = N))

d.annot_neg_ctrl_only <- d.annot_neg_ctrl_only %>%
  mutate(pgRNA_target = paste("FAKE_GENE", fake_gene_nums, sep="_"))

## get target info back into DF 
d.annot_ko <- d.annot %>%
  filter(target_type != "ctrl_ctrl") %>%
  mutate(pgRNA_target = case_when(
    target_type == "gene_gene" ~ paste(gene1_symbol, gene2_symbol, sep = "_"),
    target_type == "gene_ctrl" ~ paste(gene1_symbol, "ctrl", sep = "_"),
    target_type == "ctrl_gene" ~ paste("ctrl", gene2_symbol, sep = "_")
  ))

## - add ctrl_ctrl pgRNAs back into the df
d.annot_ctrl <- d.annot %>%
  filter(target_type == "ctrl_ctrl") %>%
  left_join(d.annot_neg_ctrl_only, by = "pgRNA_id") 

d.annot <- bind_rows(d.annot_ko, d.annot_ctrl)

#Alice had to fix
## remove gene symbols and info for control sgRNAs in single-targeting pgRNAs

d.annot <- d.annot %>%
  mutate_at(vars(starts_with("gene1")), ~replace(., target_type == "ctrl_gene", NA)) %>%
  mutate_at(vars(starts_with("gene2")), ~replace(., target_type == "gene_ctrl", NA))

## rearrange columns in d.annot
d.annot %>% dplyr::select(any_of(c("pgRNA_id", "paralog_pair", "target_type", "pgRNA_target", "gRNA1_seq", "gRNA2_seq", "norm_ctrl_flag", "paralog_pair_id")), starts_with("gene1"), starts_with("gene2"))

save_tbl(d.annot)

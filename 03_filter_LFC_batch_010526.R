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
cran_packages <- c("tidyverse", "tidylog", "kableExtra", "pheatmap", "RColorBrewer", "corrr")
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
library(RColorBrewer) # for heatmap colors
library(kableExtra) # for formatting kables
library(corrr)

args <- commandArgs(trailingOnly = TRUE)

base_dir  <- args[1]
cell_line <- args[2]

source(file.path(base_dir, "00-shared_functions_and_variables.R"))

contour_palette <- colorRampPalette(brewer.pal(n = 9, name ="Spectral"))(50)

## file paths
in_dir <- file.path(base_dir, "results", "pgRNA_counts_QC")

annot_dir <- file.path(base_dir, "results", "pgRNA_annotations")

out_dir <- file.path(base_dir, "results", "calculate_LFC")

make_out_dir(out_dir)

##functions 
make_norm_ctrl_violin_plot <- function(df, y_var, y_lab){
  
  plot <- ggplot(df, aes(x = norm_ctrl_flag, y = get(y_var), fill = norm_ctrl_flag)) +
    geom_hline(yintercept = 0) +
    geom_violin() +
    geom_boxplot(outlier.shape = NA, coef = 0, width = 0.1) +
    labs(x = "pgRNA_category", y = y_lab) +
    plot_options +
    plot_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          aspect.ratio = wide_ar) 
  
  return(plot)
  
}

get_rep_cor_dfs <- function(df, group_var, plot_var){
  
  d.reps <- df %>%
    ungroup() %>%
    dplyr::select(rep) %>%
    distinct()
  
  d.var_rep_cor <- df %>%
    ungroup() %>%
    dplyr::select(group_var, rep, plot_var) %>%
    pivot_wider(names_from = rep,
                values_from = plot_var) 
  
  print(dim(d.var_rep_cor))    # Shows number of rows and columns
  print(head(d.var_rep_cor))    # Shows first several rows
  
  d.cor <- d.var_rep_cor %>%
    dplyr::select(-group_var) %>%
    corrr::correlate() %>%
    shave() %>%
    stretch() %>%
    filter(!is.na(r)) %>%
    unite(c(x, y), col = "comparison", sep = "_vs_", remove = FALSE) %>%
    rename("sample1" = x, "sample2" = y)
  
  print(nrow(d.cor))
  print(nrow(d.var_rep_cor))
  print(sum(complete.cases(d.var_rep_cor[ , -which(names(d.var_rep_cor) == group_var)])))
  
  comparisons <- d.cor %>% pull(comparison)
  
  results <- lapply(comparisons, function(i){
    
    sample1 <- d.cor %>% filter(comparison == i) %>% pull(sample1)
    sample2 <- d.cor %>% filter(comparison == i) %>% pull(sample2)
    
    d.comparison_var <- d.var_rep_cor %>%
      dplyr::select(sample1, sample2) %>%
      rename("sample1" = sample1, "sample2" = sample2) %>%
      mutate("comparison" = i)
  })
  
  d.var_rep_cor_plot <- bind_rows(results)
  
  out_list <- list("plot_df" = d.var_rep_cor_plot, "cor_df" = d.cor)
  return(out_list)
  
}


make_rep_cor_plot <- function(plot_df, cor_df, axis_label){
  
  plot <- ggplot(plot_df, aes(x = sample1, y = sample2)) + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_point(size = 1) +
    geom_density2d(aes(color = ..level..), alpha = 0.7) +
    geom_text(data = cor_df, 
              mapping = aes(x = -Inf, y = Inf, label = paste0("r=", round(r, 3))), 
              hjust = -0.25, vjust = 1.75, size = 3.5) +
    scale_colour_gradientn(colors = rev(contour_palette)) +
    labs(x = paste("sample1", axis_label, sep = "_"),
         y = paste("sample2", axis_label, sep = "_"), 
         color = "density") +
    plot_options +
    plot_theme +
    theme(aspect.ratio = square_ar) +
    facet_wrap(~comparison)
  
  return(plot)
  
}

## counts file, flagged with whether or not to keep each pgRNA
#Alice edit to tsv
d.counts_flagged <- read_tsv(file.path(in_dir, "tables", "tsv", paste0(cell_line, "_counts_cpm_flag_long.txt")))

## annotations file
d.annot <- read_rds(file.path(annot_dir, "tables", "rds", paste0("d.", cell_line, "_annot")))


#New code by Alice to parse days and reps from current naming convention
d.counts_flagged_timepoint <-
  d.counts_flagged %>%
  separate(sample, into = c("cell", "rest"), sep = "Day", remove = FALSE) %>%
  separate(rest, into = c("day_str", "rep"), sep = "_") %>%
  mutate(
    day = readr::parse_number(day_str),    # "00" -> 0
    # rep is already "RepA" from the second separate
  ) %>%
  filter(cell !="pDNA")



d.counts_flagged_timepoint_named <- d.counts_flagged_timepoint %>%
  mutate(timepoint = case_when(
    day == min(day) ~ "T0",
    day != 0 ~ "late" )
  )

d.counts_flagged_timepoint_named

## calculate LFC for late vs T0 samples, using one T0 replicate and different late replicates
d.T0 <- d.counts_flagged_timepoint_named %>%
  filter(timepoint == "T0") %>%
  group_by(id) %>%
  summarize(T0_log2_cpm = mean(log2_cpm), .groups = "drop")

d.late_lfc <- d.counts_flagged_timepoint_named %>%
  dplyr::filter(timepoint == "late") %>%
  dplyr::select(id, sample, day, rep, count, late_log2_cpm = log2_cpm, low_T0_cpm, zero_count) %>%
  dplyr::left_join(d.T0, by = "id") %>%
  dplyr::mutate(lfc_late_vs_T0 = late_log2_cpm - T0_log2_cpm)  

d.T0
d.late_lfc

##I am up to here at the end of 12/30 - need to update pgMI mapping to double targeting/single targeting - see lfc.annot file
##Alice updated 12/18/15 to change filter to only T0 filter

d.lfc_annot <- d.late_lfc %>%
  left_join(d.annot, by = c("id" = "pgRNA_id")) %>%
  rename("pgRNA_id" = id)

d.removed_pgRNAs <- d.lfc_annot %>%
  dplyr::select(pgRNA_id, paralog_pair, low_T0_cpm) %>%
  distinct(pgRNA_id, .keep_all = TRUE) %>%
  filter(low_T0_cpm == TRUE) 

d.removed_pgRNAs
save_tbl(d.removed_pgRNAs)

##filter out pgRNAs to remove

d.lfc_annot <- d.lfc_annot %>%
  filter(low_T0_cpm == FALSE)

d.lfc_annot

save_tbl(d.lfc_annot)

unadjusted_norm_ctrl_violin_plot <- make_norm_ctrl_violin_plot(d.lfc_annot, "lfc_late_vs_T0", "raw_LFC")

unadjusted_norm_ctrl_violin_plot +
  facet_wrap(~rep)

save_plot(unadjusted_norm_ctrl_violin_plot)



d.lfc_annot %>%
  filter(norm_ctrl_flag == "negative_control" | norm_ctrl_flag == "positive_control") %>%
  group_by(rep, norm_ctrl_flag) %>%
  summarize(mean = mean(lfc_late_vs_T0),
            sd = sd(lfc_late_vs_T0)) %>%
  pivot_wider(names_from = norm_ctrl_flag, 
              values_from = c(mean, sd),
              names_glue = "{norm_ctrl_flag}_{.value}") %>%
  mutate(ssmd = (positive_control_mean - negative_control_mean) /
           (sqrt(positive_control_sd^2 + negative_control_sd^2)))

##calculate SSMD

d.lfc_annot %>%
  ungroup() %>%
  filter(norm_ctrl_flag == "negative_control" | norm_ctrl_flag == "positive_control") %>%
  mutate(keep_flag = case_when(
    target_type == "ctrl_ctrl" ~ TRUE,
    target_type == "gene_ctrl" & gene1_expressed_flag == TRUE ~ TRUE,
    target_type == "ctrl_gene" & gene2_expressed_flag == TRUE ~ TRUE,
    TRUE ~ FALSE)) %>% ## if these conditions are not met, set to F
  filter(keep_flag == TRUE) %>%
  group_by(rep, norm_ctrl_flag) %>%
  summarize(mean = mean(lfc_late_vs_T0),
            sd = sd(lfc_late_vs_T0)) %>%
  pivot_wider(names_from = norm_ctrl_flag, 
              values_from = c(mean, sd),
              names_glue = "{norm_ctrl_flag}_{.value}") %>%
  mutate(ssmd = (positive_control_mean - negative_control_mean) /
           (sqrt(positive_control_sd^2 + negative_control_sd^2)))

d.ssmd <- d.lfc_annot %>%
  filter(norm_ctrl_flag == "negative_control" | norm_ctrl_flag == "positive_control") %>%
  group_by(rep, norm_ctrl_flag) %>%
  summarize(mean = mean(lfc_late_vs_T0),
            var = var(lfc_late_vs_T0)) %>%
  pivot_wider(names_from = norm_ctrl_flag, 
              values_from = c(mean, var),
              names_glue = "{norm_ctrl_flag}_{.value}") %>%
  mutate(ssmd = (positive_control_mean - negative_control_mean) /
           (sqrt(positive_control_var + negative_control_var)))
print(d.ssmd)

save_tbl(d.ssmd)

##Alice add NNMD
d.nnmd <- d.lfc_annot %>%
  filter(norm_ctrl_flag == "negative_control" | norm_ctrl_flag == "positive_control") %>%
  group_by(rep, norm_ctrl_flag) %>%
  summarize(median = median(lfc_late_vs_T0),
            mad = mad(lfc_late_vs_T0)) %>%
  pivot_wider(names_from = norm_ctrl_flag, 
              values_from = c(median, mad),
              names_glue = "{norm_ctrl_flag}_{.value}") %>%
  mutate(nnmd = (positive_control_median - negative_control_median) /
           negative_control_mad)
print(d.nnmd)
save_tbl(d.nnmd)

dim(d.lfc_annot)

d.lfc_annot_no_dups <- d.lfc_annot %>%
  distinct(pgRNA_id, rep, .keep_all = TRUE)

dim(d.lfc_annot_no_dups)

results <- get_rep_cor_dfs(d.lfc_annot_no_dups, "pgRNA_id", "lfc_late_vs_T0")

d.unadj_lfc_rep_cor_plot <- results[[1]]
d.unadj_lfc_rep_cor_summary <- results[[2]]


print_kbl(d.unadj_lfc_rep_cor_summary)

## make replicate correlation plot
make_rep_cor_plot(d.unadj_lfc_rep_cor_plot, d.unadj_lfc_rep_cor_summary, "LFC")
save_plot(make_rep_cor_plot(d.unadj_lfc_rep_cor_plot, d.unadj_lfc_rep_cor_summary, "LFC"))

##adjust LFCs

d.control_group_medians <- d.lfc_annot %>%
  group_by(rep, norm_ctrl_flag) %>%
  filter(norm_ctrl_flag == "negative_control" | norm_ctrl_flag == "positive_control") %>%
  summarize(median_lfc = median(lfc_late_vs_T0))
print_kbl(d.control_group_medians)

ctrl_medians <- d.lfc_annot %>%
  filter(norm_ctrl_flag == "negative_control") %>%
  group_by(rep) %>%
  summarise(neg_ctrl_median = median(lfc_late_vs_T0, na.rm = TRUE))

pos_ctrl_medians <- d.lfc_annot %>%
  filter(norm_ctrl_flag == "positive_control") %>%
  group_by(rep) %>%
  summarise(pos_ctrl_median = median(lfc_late_vs_T0, na.rm = TRUE))

d.lfc_annot_adj <- d.lfc_annot %>%
  left_join(ctrl_medians, by = "rep") %>%
  left_join(pos_ctrl_medians, by = "rep") %>%
  mutate(
    lfc_adj1 = lfc_late_vs_T0 - neg_ctrl_median,
    lfc_adj2 = lfc_adj1 / (neg_ctrl_median - pos_ctrl_median)
  )


CS_pgRNA_violin_plot <- make_norm_ctrl_violin_plot(d.lfc_annot_adj, "lfc_adj2", "CRISPR score")

CS_pgRNA_violin_plot +
    facet_wrap(~rep)

save_plot(CS_pgRNA_violin_plot)
save_plot(CS_pgRNA_violin_plot + facet_wrap((~rep)))

d.lfc_annot_adj_single <- d.lfc_annot_adj %>%
  filter(target_type == "gene_ctrl" | target_type == "ctrl_gene") %>%
  ## make a flag variable to indicate which pgRNAs are targeting unexpressed
  ## single targets
  mutate(unexpressed_ctrl_flag = case_when(
    target_type == "gene_ctrl" & gene1_expressed_flag == FALSE ~ TRUE,
    target_type == "ctrl_gene" & gene2_expressed_flag == FALSE ~ TRUE,
    TRUE ~ FALSE 
  )) 
## %>% ##Alice commented out so we could still flag/summarize but not adjust
##  group_by(rep) %>%
##  mutate(lfc_adj3 = lfc_adj2 - median(lfc_adj2[unexpressed_ctrl_flag == TRUE]))

d.lfc_annot_adj_single_summary <- d.lfc_annot_adj_single %>%
  group_by(rep, unexpressed_ctrl_flag) %>%
  summarize(median = median(lfc_adj2)) ##switched to adj2
print_kbl(d.lfc_annot_adj_single_summary)
##this is now summarizing the median adj lfc depending on expression of single targeting

##Alice adjusted to keep showing adj2 (CRISPR score) instead of adj3 (CS adjusted by Exp)

CS_by_expression_singletargeting <- d.lfc_annot_adj_single %>%
  filter(!is.na(gene1_log2_tpm) | !is.na(gene2_log2_tpm)) %>%
  ggplot(aes(x = lfc_adj2, fill = unexpressed_ctrl_flag)) +
  geom_density(alpha = 0.7) +
  geom_vline(data = d.lfc_annot_adj_single_summary, 
             aes(xintercept = median, color = unexpressed_ctrl_flag),
             linetype = "dashed") +
  scale_fill_discrete(name = "expression_group", 
                      limits = c(FALSE, TRUE),
                      labels = c("expressed", "unexpressed")) +
  scale_color_discrete(name = "expression_group", 
                       limits = c(FALSE, TRUE),
                       labels = c("expressed", "unexpressed")) +
  labs(x = "CRISPR score") +
  plot_options +
  plot_theme +
  theme(aspect.ratio = wide_ar) +
  facet_wrap(~rep)

CS_by_expression_singletargeting
save_plot(CS_by_expression_singletargeting)

d.lfc_annot_adj_double <- d.lfc_annot_adj %>%
  filter(target_type == "gene_gene") %>%
  ## make a flag variable to indicate which pgRNAs are targeting double
  ## unexpressed targets
  mutate(unexpressed_ctrl_flag = case_when(
    gene1_expressed_flag == FALSE & gene2_expressed_flag == FALSE ~ TRUE,
    TRUE ~ FALSE)) 
## %>% ##Alice commented out again to allow summary but not adjustment
##  group_by(rep) %>%
##  mutate(lfc_adj3 = lfc_adj2 - median(lfc_adj2[unexpressed_ctrl_flag == TRUE]))

d.lfc_annot_adj_double_summary <- d.lfc_annot_adj_double %>%
  group_by(rep, unexpressed_ctrl_flag) %>%
  summarize(median = median(lfc_adj2))
print_kbl(d.lfc_annot_adj_double_summary)

d.lfc_annot_adj_double_plot <-  d.lfc_annot_adj_double %>% 
  filter(!is.na(gene1_log2_tpm) & !is.na(gene2_log2_tpm)) %>%
  mutate(n_genes_expressed = case_when(
    gene1_expressed_flag == FALSE & gene2_expressed_flag == FALSE ~ "0",
    gene1_expressed_flag == TRUE & gene2_expressed_flag == FALSE ~ "1",
    gene1_expressed_flag == FALSE & gene2_expressed_flag == TRUE ~ "1",
    gene1_expressed_flag == TRUE & gene2_expressed_flag == TRUE ~ "2")) 

d.lfc_annot_adj_double_plot_summary <- d.lfc_annot_adj_double_plot %>%
  group_by(rep, n_genes_expressed) %>%
  summarize(median = median(lfc_adj2))


CS_by_expression_doubletargeting <- ggplot(d.lfc_annot_adj_double_plot, aes(x = lfc_adj2, fill = n_genes_expressed)) +
  geom_density(alpha = 0.7) +
  geom_vline(data = d.lfc_annot_adj_double_plot_summary,
             aes(xintercept = median, color = n_genes_expressed),
             linetype = "dashed") +
  labs(x = "CRISPR score") +
  plot_options +
  plot_theme +
  theme(aspect.ratio = wide_ar) +
  facet_wrap(~rep)

CS_by_expression_doubletargeting
save_plot(CS_by_expression_doubletargeting)

##Left of here 12/31 - need to figure out what rejoining is required if any
### ntc_ntc
d.lfc_annot_adj_control <- d.lfc_annot_adj %>%
  filter(target_type == "ctrl_ctrl") %>%
  mutate(lfc_adj3 = lfc_adj2)

### single targeting
# colnames(d.lfc_annot_adj_single)

d.lfc_annot_adj_single <- d.lfc_annot_adj_single %>%
  dplyr::select(-unexpressed_ctrl_flag)

### double targeting
d.lfc_annot_adj_double <- d.lfc_annot_adj_double %>%
  dplyr::select(-unexpressed_ctrl_flag)

## bind rows
d.lfc_annot_adj_pgRNA <- bind_rows(d.lfc_annot_adj_double, d.lfc_annot_adj_single, d.lfc_annot_adj_control)

## rename final adjusted column to CRISPR_score
d.lfc_annot_adj_pgRNA <- d.lfc_annot_adj_pgRNA %>%
  rename(CRISPR_score = lfc_adj2) %>% ##Alice updated back to use adj2 
  dplyr::select(-lfc_adj1)

save_tbl(d.lfc_annot_adj_pgRNA)

##target level values
d.lfc_annot_adj_target <- d.lfc_annot_adj_pgRNA %>%
  group_by(rep, pgRNA_target) %>%
  mutate(target_mean_CS = mean(CRISPR_score),
         target_median_CS = median(CRISPR_score)) %>%
  distinct(pgRNA_target, .keep_all = TRUE) %>%
  dplyr::select(-c(pgRNA_id, CRISPR_score, contains("seq")))

#save rep level target-level CRISPR scores  
save_tbl(d.lfc_annot_adj_target)

results <- get_rep_cor_dfs(d.lfc_annot_adj_target, "pgRNA_target", "target_mean_CS")

d.adj_mean_lfc_rep_cor_plot <- results[[1]]
d.adj_mean_lfc_rep_cor_summary <- results[[2]]

print_kbl(d.adj_mean_lfc_rep_cor_summary)

## make replicate correlation plot
rep_cor_plot <- make_rep_cor_plot(d.adj_mean_lfc_rep_cor_plot, d.adj_mean_lfc_rep_cor_summary, "target_mean_CS")

rep_cor_plot
save_plot(rep_cor_plot)

## target-level violin plot
CS_target_violin_plot <- make_norm_ctrl_violin_plot(d.lfc_annot_adj_target, "target_mean_CS", "target_mean_CRISPR_score")

CS_target_violin_plot +
    facet_wrap(~rep)

save_plot(CS_target_violin_plot)
save_plot(CS_target_violin_plot + facet_wrap(~rep))

## mean across reps violin plot
d.lfc_annot_adj_target_rep_mean <- d.lfc_annot_adj_target %>%
  group_by(pgRNA_target) %>%
  mutate(rep_target_mean_CS = mean(target_mean_CS)) %>%
  distinct(pgRNA_target, .keep_all = TRUE) %>%
  dplyr::select(-c(rep, target_mean_CS, target_median_CS))


rep_mean_CS_target_violin_plot <- make_norm_ctrl_violin_plot(d.lfc_annot_adj_target_rep_mean, "rep_target_mean_CS", "mean_targetCS_across_reps")

rep_mean_CS_target_violin_plot
save_plot(rep_mean_CS_target_violin_plot)

## get (& plot) mean across retained reps
d.lfc_annot_adj_pgRNA_rep_mean <- d.lfc_annot_adj_pgRNA %>%
  group_by(pgRNA_id) %>%
  mutate(rep_mean_CS = mean(CRISPR_score)) %>%
  distinct(pgRNA_target, .keep_all = TRUE) %>%
  dplyr::select(-c(rep, CRISPR_score))


rep_mean_CS_pgRNA_violin_plot <- make_norm_ctrl_violin_plot(d.lfc_annot_adj_pgRNA_rep_mean, "rep_mean_CS", "mean_pgRNA_CS")

rep_mean_CS_pgRNA_violin_plot
save_plot(rep_mean_CS_pgRNA_violin_plot)




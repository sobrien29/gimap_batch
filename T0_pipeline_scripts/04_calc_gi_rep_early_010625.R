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
cran_packages <- c("tidyverse", "tidylog", "kableExtra", "pheatmap", "RColorBrewer", "ggrepel")
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

args <- commandArgs(trailingOnly = TRUE)

base_dir  <- args[1]
cell_line <- args[2]


source(file.path(base_dir, "00-shared_functions_and_variables.R"))

contour_palette <- colorRampPalette(brewer.pal(n = 9, name ="Spectral"))(50)

## file paths
in_dir <- file.path(base_dir, "results", "calculate_LFC")

out_dir <- file.path(base_dir, "results", "calculate_GI_scores", "pre-rep_averaging")

make_out_dir(out_dir)

## counts file, flagged with whether or not to keep each pgRNA
d.lfc_pgRNA <- read_rds(file.path(in_dir, "tables", "rds", paste0("d.", cell_line, "_lfc_annot_adj_pgRNA")))

##lets average the CRISPR scores of all three replicates

d.lfc_pgRNA_avg <- d.lfc_pgRNA %>%
  dplyr::group_by(pgRNA_id) %>%
  dplyr::summarize(
    CRISPR_score_avg = mean(CRISPR_score, na.rm = TRUE),
    
    # keep pgRNA-level metadata (must be constant within pgRNA_id)
    paralog_pair = dplyr::first(paralog_pair),
    gRNA1_seq = dplyr::first(gRNA1_seq),
    gRNA2_seq = dplyr::first(gRNA2_seq),
    target_type = dplyr::first(target_type),
    pgRNA_target = dplyr::first(pgRNA_target),
    gene1_symbol = dplyr::first(gene1_symbol),
    gene2_symbol = dplyr::first(gene2_symbol),
  
    
    .groups = "drop"
  )

## get just double-targeting pgRNAs
d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA_avg %>%
  ungroup() %>%
  filter(target_type == "gene_gene") %>%
  dplyr::select(pgRNA_id, paralog_pair, CRISPR_score_avg, target_type, pgRNA_target,
                gRNA1_seq, gRNA2_seq, gene1_symbol, gene2_symbol)



## calculate mean CRISPR score of single-targeting pgRNAs containing the same targeting
## sgRNA sequence but different control sgRNA sequences
d.mean_single_target_CS <- d.lfc_pgRNA_avg %>%
  ungroup() %>%
  filter(target_type == "gene_ctrl" | target_type == "ctrl_gene") %>%
  mutate(targeting_gRNA_seq = case_when(
    target_type == "gene_ctrl" ~ gRNA1_seq,
    target_type == "ctrl_gene" ~ gRNA2_seq
  )) %>%
  group_by( paralog_pair, targeting_gRNA_seq) %>%
  mutate(mean_single_target_CS = mean(CRISPR_score_avg)) %>%
  dplyr::select(paralog_pair, targeting_gRNA_seq, mean_single_target_CS) %>%
  distinct(paralog_pair, targeting_gRNA_seq, .keep_all = TRUE)



## join single-target CRISPR scores with double-targeting pgRNA df based on targeting
## sgRNA sequences
d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA_double_targeting %>%
  rename(double_target_CS = CRISPR_score_avg) %>%
  left_join(d.mean_single_target_CS, by = c("paralog_pair", "gRNA1_seq" = "targeting_gRNA_seq")) %>%
  rename(mean_gRNA1_single_target_CS = mean_single_target_CS) %>%
  left_join(d.mean_single_target_CS, by = c("paralog_pair", "gRNA2_seq" = "targeting_gRNA_seq")) %>%
  rename(mean_gRNA2_single_target_CS = mean_single_target_CS)

##lets remove NAs here for pgRNAs that were removed from the data
d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA_double_targeting %>% 
  filter(!is.na(mean_gRNA1_single_target_CS)) %>% 
  filter(!is.na(mean_gRNA2_single_target_CS))

## calculate expected double-targeting GI score by summing the two mean single-targeting
## CRISPR scores for that paralog pair
d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA_double_targeting %>% mutate(expected_CS = mean_gRNA1_single_target_CS + mean_gRNA2_single_target_CS)

## get just single-targeting pgRNAs
d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_avg %>%
  ungroup() %>%
  filter(target_type == "gene_ctrl" | target_type == "ctrl_gene") %>%
  dplyr::select(pgRNA_id,  paralog_pair, CRISPR_score_avg, target_type, pgRNA_target,
                gRNA1_seq, gRNA2_seq, gene1_symbol, gene2_symbol) %>%
  mutate(targeting_gRNA_seq = case_when(
    target_type == "gene_ctrl" ~ gRNA1_seq,
    target_type == "ctrl_gene" ~ gRNA2_seq
  )) %>%
  mutate(control_gRNA_seq = case_when(
    target_type == "gene_ctrl" ~ gRNA2_seq,
    target_type == "ctrl_gene" ~ gRNA1_seq
  ))

## get targeting gRNAs to add back to DF to calculate expected GI scores
d.other_single_targeting_CS <- d.lfc_pgRNA_single_targeting %>%
  dplyr::select(paralog_pair, targeting_gRNA_seq, control_gRNA_seq, CRISPR_score_avg) %>%
  rename(other_single_target_CS = CRISPR_score_avg, other_control_seq = control_gRNA_seq)

## add back "other" single-targeting pgRNA CRISPR scores into the main DF, 
## get rid of duplicates (same target seq and same control seq)
d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  left_join(d.other_single_targeting_CS, by = c("paralog_pair", "targeting_gRNA_seq")) %>%
  mutate(same_control_seq = ifelse(control_gRNA_seq == other_control_seq, TRUE, FALSE)) %>%
  filter(same_control_seq == FALSE) %>%
  dplyr::select(-same_control_seq)

## add column for double crispr controls set to 0 and generate expected CS
d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  mutate(mean_double_control_CS = 0) %>% mutate(expected_CS = other_single_target_CS + mean_double_control_CS) %>%
  rename(single_target_CS = CRISPR_score_avg)

#### Linear model
d.lfc_pgRNA_single_targeting

d.lfc_pgRNA_single_targeting_mean <- d.lfc_pgRNA_single_targeting %>%
  group_by(pgRNA_target) %>%
  summarize(mean_expected_CS = mean(expected_CS),
            mean_observed_CS = mean(single_target_CS)) %>%
  ungroup()
d.lfc_pgRNA_single_targeting_mean

## fit linear model to target-level mean single-targeting pgRNA expected vs. 
## observed values and extract slope and intercept values
d.lfc_pgRNA_single_targeting_mean_lm_summary <- d.lfc_pgRNA_single_targeting_mean %>%
  group_modify(~ broom::tidy(lm(mean_observed_CS ~ mean_expected_CS, data = .x))) %>%
  dplyr::ungroup() %>%
  dplyr::select(term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(intercept = "(Intercept)", slope = mean_expected_CS)
print_kbl(d.lfc_pgRNA_single_targeting_mean_lm_summary)

## plot
d.lfc_pgRNA_single_targeting_mean.plot <- d.lfc_pgRNA_single_targeting_mean %>%
  ggplot(aes(x = mean_expected_CS, y = mean_observed_CS)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "gray55") +
  labs(x = "expected_CRISPR_score", y = "observed_CRISPR_score") +
  plot_options +
  plot_theme +
  theme(aspect.ratio = square_ar) 
d.lfc_pgRNA_single_targeting_mean.plot
save_plot(d.lfc_pgRNA_single_targeting_mean.plot)

d.lfc_single_targeting_expected_meanobserved <- ggplot(data = d.lfc_pgRNA_single_targeting, 
                                                       aes(x = expected_CS, y = single_target_CS)) +
  geom_point() + 
  geom_density2d(aes(color = ..level..), alpha = 0.9) +
  geom_abline(data = d.lfc_pgRNA_single_targeting_mean_lm_summary,
              aes(slope = slope, intercept = intercept), 
              color = "gray55", size = 1) +
  scale_colour_gradientn(colors = rev(contour_palette)) +
  labs(x = "expected_CRISPR_score", 
       y = "observed_CRISPR_score",
       color = "density") +
  plot_options +
  plot_theme +
  theme(aspect.ratio = square_ar) 
d.lfc_single_targeting_expected_meanobserved
save_plot(d.lfc_single_targeting_expected_meanobserved)

d.lfc_single_targeting_expected_meanobserved_count <- ggplot(data = d.lfc_pgRNA_single_targeting, 
                                                             aes(x = expected_CS, y = single_target_CS)) +
  geom_point() + 
  geom_density2d(aes(color = ..level..), alpha = 0.9, contour_var = "count") +
  geom_abline(data = d.lfc_pgRNA_single_targeting_mean_lm_summary,
              aes(slope = slope, intercept = intercept), 
              color = "gray55", size = 1) +
  scale_colour_gradientn(colors = rev(contour_palette)) +
  labs(x = "expected_CRISPR_score", 
       y = "observed_CRISPR_score",
       color = "count") +
  plot_options +
  plot_theme +
  theme(aspect.ratio = square_ar)
save_plot(d.lfc_single_targeting_expected_meanobserved_count)

### Single-target GI scores
d.lfc_pgRNA_single_targeting_GI <- d.lfc_pgRNA_single_targeting %>% 
  mutate(slope = d.lfc_pgRNA_single_targeting_mean_lm_summary$slope[1],
         intercept = d.lfc_pgRNA_single_targeting_mean_lm_summary$intercept[1]) %>%
  mutate(GI_score = single_target_CS - (intercept + slope * expected_CS)) 

# raw_GI = target_mean_CS - (control_intercept + control_slope * target_mean_expected_C
d.lfc_pgRNA_double_targeting_GI <- d.lfc_pgRNA_double_targeting %>%
  mutate(slope = d.lfc_pgRNA_single_targeting_mean_lm_summary$slope[1],
         intercept = d.lfc_pgRNA_single_targeting_mean_lm_summary$intercept[1]) %>% 
  mutate(GI_score = double_target_CS - (intercept + slope * expected_CS))

## reformat to match
d.GI_scores_pgRNA_double <- d.lfc_pgRNA_double_targeting_GI %>%
  rename(observed_CS = double_target_CS) %>%
  dplyr::select(-c(mean_gRNA1_single_target_CS, mean_gRNA2_single_target_CS)) %>%
  mutate(broad_target_type = "double_targeting") %>%
  ungroup()

d.GI_scores_pgRNA_single <- d.lfc_pgRNA_single_targeting_GI %>%
  rename(observed_CS = single_target_CS) %>%
  dplyr::select(-c(other_single_target_CS, mean_double_control_CS,
                   targeting_gRNA_seq, control_gRNA_seq, other_control_seq)) %>%
  mutate(broad_target_type = "single_targeting") %>%
  ungroup()

### Calculate p-values
## t-test (parametric)
## get a vector of GI scores for all single-targeting ("control") pgRNAs
single_GI_scores <- d.GI_scores_pgRNA_single %>%
  pull(GI_score)

## for each paralog pair, test double-targeting GI scores
d.double_GI_scores <- d.GI_scores_pgRNA_double %>%
  dplyr::group_by(paralog_pair) %>%
  dplyr::mutate(
    p_val = t.test(
      x = single_GI_scores,
      y = GI_score,
      paired = FALSE
    )$p.value
  )

## adjust for multiple testing using Benjamini–Hochberg
d.p_val <- d.double_GI_scores %>%
  dplyr::select(paralog_pair, p_val) %>%
  dplyr::arrange(p_val) %>%
  dplyr::distinct(p_val, .keep_all = TRUE)

fdr_vals <- p.adjust(d.p_val$p_val, method = "BH")

d.fdr <- tibble(fdr = fdr_vals) %>%
  dplyr::bind_cols(d.p_val) %>%
  dplyr::select(paralog_pair, fdr)

## add FDR values back into the double-targeting DF
results <- d.double_GI_scores %>%
  dplyr::left_join(d.fdr, by = "paralog_pair")



d.GI_scores_pgRNA_double <- bind_rows(results)

### Bind double- and single-targeting dfs
d.stats <- d.GI_scores_pgRNA_double %>%
  dplyr::select(paralog_pair, p_val, fdr) %>%
  distinct(paralog_pair, .keep_all = TRUE) 

## add p-val and fdr to single-targeting
d.GI_scores_pgRNA_single <- d.GI_scores_pgRNA_single %>%
  left_join(d.stats, by = "paralog_pair")

d.GI_scores_pgRNA <- bind_rows(d.GI_scores_pgRNA_double, d.GI_scores_pgRNA_single)

d.GI_scores_target <- d.GI_scores_pgRNA %>%
  ungroup() %>%
  group_by(pgRNA_target) %>%
  mutate(mean_observed_CS = mean(observed_CS),
         mean_expected_CS = mean(expected_CS),
         mean_GI_score = mean(GI_score)) %>%
  distinct(pgRNA_target, .keep_all = TRUE) %>%
  dplyr::select(-c(pgRNA_id, contains("seq"), observed_CS, expected_CS, GI_score)) %>%
  ungroup()

d.GI_scores_target %>%
  summarize(n = n())

d.GI_scores_target_plot <- ggplot(d.GI_scores_target, aes(x = mean_expected_CS, y = mean_observed_CS, color = broad_target_type)) +
  geom_point() + 
  geom_abline(data = d.lfc_pgRNA_single_targeting_mean_lm_summary,
              aes(slope = slope, intercept = intercept),
              color = "gray55", size = 1) +
  labs(x = "expected_CRISPR_score", y = "observed_CRISPR_score") +
  plot_options +
  plot_theme +
  theme(aspect.ratio = square_ar)
save_plot(d.GI_scores_target_plot)


## rank scatter plot 
quantiles <- quantile(d.GI_scores_target$mean_GI_score, probs = c(0.15, 0.85), na.rm = TRUE)

rank_scatter <- ggplot <- d.GI_scores_target %>%
  dplyr::ungroup() %>%
  filter(target_type == "gene_gene") %>%
  mutate(Rank = percent_rank(mean_GI_score)) %>%
  ggplot(aes(
    x = Rank,
    y = mean_GI_score
  )) +
  geom_point(size = 1, alpha = 0.7) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  ylab("GI score") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  geom_text_repel(
    data = d.GI_scores_target %>%
      dplyr::ungroup() %>%
      dplyr::filter(target_type == "gene_gene") %>%
      dplyr::mutate(Rank = percent_rank(mean_GI_score)) %>%
      dplyr::filter(
        mean_GI_score <= quantiles[1] |
          mean_GI_score >= quantiles[2]
      ),
    aes(label = paralog_pair),
    max.overlaps = 10, 
    size=3
  )

rank_scatter
save_plot(rank_scatter)
##volcano plot

volcano_plot <- gplot <- d.GI_scores_target %>%
  filter(target_type == "gene_gene") %>% # get only double targeting
  mutate(
    logfdr = -log10(fdr),
    pointColor = case_when(logfdr < 1 ~ "darkgrey",
                           ((mean_GI_score < -0.5) & (logfdr > 1)) ~ "dodgerblue3",
                           ((mean_GI_score > 0.5) & (logfdr > 1)) ~ "darkred",
                           .default = "black"
    )
  ) %>%
  ggplot(aes(
    x = mean_GI_score,
    y = logfdr,
    color = pointColor
  )) +
  geom_point(size = 1, alpha = 0.7) +
  theme_classic() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  theme(legend.position = "none") +
  scale_color_manual(values = c(
    "darkgrey" = "darkgrey",
    "dodgerblue3" = "dodgerblue3",
    "darkred" = "darkred",
    "black" = "black"
  )) +
  ylab("-log10(FDR)") +
  xlab("Mean GI score") +
  geom_text_repel(
    data = function(df) {
      df %>% filter(pointColor %in% c("dodgerblue3", "darkred"))
    },
    aes(label = paralog_pair),
    size = 3,
    max.overlaps = 15,
    box.padding = 0.4,
    point.padding = 0.3,
    show.legend = FALSE
  )

volcano_plot
save_plot(volcano_plot)

## Save output
save_tbl(d.GI_scores_target)
save_tbl(d.GI_scores_pgRNA)

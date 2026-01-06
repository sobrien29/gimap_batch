library(tidyverse)
library(tidylog)
library(RColorBrewer) # for heatmap colors
library(kableExtra) # for formatting kables

source(file.path(base_dir, "00-shared_functions_and_variables.R"))

contour_palette <- colorRampPalette(brewer.pal(n = 9, name ="Spectral"))(50)

## file paths
in_dir <- file.path(base_dir, "results", "calculate_LFC")

out_dir <- file.path(base_dir, "results", "calculate_GI_scores")

make_out_dir(out_dir)

## counts file, flagged with whether or not to keep each pgRNA
d.lfc_pgRNA <- read_rds(file.path(in_dir, "tables", "rds", paste0("d.", cell_line, "_lfc_annot_adj_pgRNA")))


## get just double-targeting pgRNAs
d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "gene_gene") %>%
  dplyr::select(pgRNA_id, rep, paralog_pair, CRISPR_score, target_type, pgRNA_target,
                gRNA1_seq, gRNA2_seq, gene1_symbol, gene2_symbol)

## confirm that getting mean single-targeting CRISPR scores is acting as expected:
d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "gene_ctrl" | target_type == "ctrl_gene") %>%
  mutate(targeting_gRNA_seq = case_when(
    target_type == "gene_ctrl" ~ gRNA1_seq,
    target_type == "ctrl_gene" ~ gRNA2_seq
  )) %>%
  group_by(rep, paralog_pair, targeting_gRNA_seq) %>%
  mutate(mean_single_target_CS = mean(CRISPR_score)) %>%
  dplyr::select(pgRNA_id:gRNA2_seq, targeting_gRNA_seq, mean_single_target_CS) %>%
  arrange(rep, paralog_pair, targeting_gRNA_seq)

## calculate mean CRISPR score of single-targeting pgRNAs containing the same targeting
## sgRNA sequence but different control sgRNA sequences
d.mean_single_target_CS <- d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "gene_ctrl" | target_type == "ctrl_gene") %>%
  mutate(targeting_gRNA_seq = case_when(
    target_type == "gene_ctrl" ~ gRNA1_seq,
    target_type == "ctrl_gene" ~ gRNA2_seq
  )) %>%
  group_by(rep, paralog_pair, targeting_gRNA_seq) %>%
  mutate(mean_single_target_CS = mean(CRISPR_score)) %>%
  dplyr::select(rep, paralog_pair, targeting_gRNA_seq, mean_single_target_CS) %>%
  distinct(rep, paralog_pair, targeting_gRNA_seq, .keep_all = TRUE)

## add mean single-targeting CRISPR scores into double-targeting DF so I can calculate
## expected GI scores by summing 

## join single-target CRISPR scores with double-targeting pgRNA df based on targeting
## sgRNA sequences
d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA_double_targeting %>%
  rename(double_target_CS = CRISPR_score) %>%
  left_join(d.mean_single_target_CS, by = c("rep", "paralog_pair", "gRNA1_seq" = "targeting_gRNA_seq")) %>%
  rename(mean_gRNA1_single_target_CS = mean_single_target_CS) %>%
  left_join(d.mean_single_target_CS, by = c("rep", "paralog_pair", "gRNA2_seq" = "targeting_gRNA_seq")) %>%
  rename(mean_gRNA2_single_target_CS = mean_single_target_CS)

## save a list of pgRNAs that are filtered out based on single targeting CS being NA
d.rm_double_targeting_pgRNAs <- d.lfc_pgRNA_double_targeting %>%
  dplyr::select(pgRNA_id, rep, paralog_pair, double_target_CS, gRNA1_seq, gRNA2_seq, 
                mean_gRNA1_single_target_CS, mean_gRNA2_single_target_CS) %>%
  filter(is.na(mean_gRNA1_single_target_CS) | is.na(mean_gRNA2_single_target_CS))
save_tbl(d.rm_double_targeting_pgRNAs)

## filter out single-targeting CS == NA rows
d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA_double_targeting %>%
  ## filter out pgRNAs where either single-targeting mean CRISPR score is NA
  filter(!is.na(mean_gRNA1_single_target_CS) & !is.na(mean_gRNA2_single_target_CS)) %>%
  ## calculate expected double-targeting GI score by summing the two mean single-targeting
  ## CRISPR scores for that paralog pair
  mutate(expected_CS = mean_gRNA1_single_target_CS + mean_gRNA2_single_target_CS)

## get just single-targeting pgRNAs
d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "gene_ctrl" | target_type == "ctrl_gene") %>%
  dplyr::select(pgRNA_id, rep, paralog_pair, CRISPR_score, target_type, pgRNA_target,
                gRNA1_seq, gRNA2_seq, gene1_symbol, gene2_symbol) %>%
  mutate(targeting_gRNA_seq = case_when(
    target_type == "gene_ctrl" ~ gRNA1_seq,
    target_type == "ctrl_gene" ~ gRNA2_seq
  )) %>%
  mutate(control_gRNA_seq = case_when(
    target_type == "gene_ctrl" ~ gRNA2_seq,
    target_type == "ctrl_gene" ~ gRNA1_seq
  ))

## calculate mean CRISPR score of double-non-targeting pgRNAs containing the same NTC
## sgRNA sequence
d.mean_double_control_CS <- d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "ctrl_ctrl") %>%
  pivot_longer(cols = c(gRNA1_seq, gRNA2_seq),
               names_to = "position", 
               values_to = "control_gRNA_seq") %>%
  group_by(rep, control_gRNA_seq) %>%
  mutate(mean_double_control_CS = 0) %>%
  dplyr::select(rep, control_gRNA_seq, mean_double_control_CS) %>%
  distinct(rep, control_gRNA_seq, .keep_all = TRUE)

## get targeting gRNAs to add back to DF to calculate expected GI scores
d.other_single_targeting_CS <- d.lfc_pgRNA_single_targeting %>%
  dplyr::select(rep, paralog_pair, targeting_gRNA_seq, control_gRNA_seq, CRISPR_score) %>%
  rename(other_single_target_CS = CRISPR_score, other_control_seq = control_gRNA_seq)

## add back "other" single-targeting pgRNA CRISPR scores into the main DF, 
## get rid of duplicates (same target seq and same control seq)
d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  left_join(d.other_single_targeting_CS, by = c("rep", "paralog_pair", "targeting_gRNA_seq")) %>%
  mutate(same_control_seq = ifelse(control_gRNA_seq == other_control_seq, TRUE, FALSE)) %>%
  filter(same_control_seq == FALSE) %>%
  dplyr::select(-same_control_seq)

## add back double control CRISPR scores too
d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  left_join(y = d.mean_double_control_CS, by = c("rep", "control_gRNA_seq"))

## save a list of pgRNAs that will be removed
d.rm_single_targeting_pgRNAs <- d.lfc_pgRNA_single_targeting %>%
  filter(is.na(other_single_target_CS) | is.na(mean_double_control_CS))
save_tbl(d.rm_single_targeting_pgRNAs)

## filter out pgRNAs whose CS = NA
d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  filter(!is.na(other_single_target_CS) & !is.na(mean_double_control_CS))

d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  mutate(expected_CS = other_single_target_CS + mean_double_control_CS) %>%
  rename(single_target_CS = CRISPR_score)

#### Linear model
d.lfc_pgRNA_single_targeting

d.lfc_pgRNA_single_targeting_mean <- d.lfc_pgRNA_single_targeting %>%
  group_by(rep, pgRNA_target) %>%
  summarize(mean_expected_CS = mean(expected_CS),
            mean_observed_CS = mean(single_target_CS)) %>%
  ungroup()
d.lfc_pgRNA_single_targeting_mean

## fit linear model to target-level mean single-targeting pgRNA expected vs. 
## observed values and extract slope and intercept values
d.lfc_pgRNA_single_targeting_mean_lm_summary <- d.lfc_pgRNA_single_targeting_mean %>%
  group_by(rep) %>%
  group_modify(~ broom::tidy(lm(mean_observed_CS ~ mean_expected_CS, data = .x))) %>%
  dplyr::ungroup() %>%
  dplyr::select(rep, term, estimate) %>%
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
  theme(aspect.ratio = square_ar) +
  facet_wrap(~rep)
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
  theme(aspect.ratio = square_ar) +
  facet_wrap(~rep)
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
  theme(aspect.ratio = square_ar) +
  facet_wrap(~rep)
save_plot(d.lfc_single_targeting_expected_meanobserved_count)


### Single-target GI scores
d.lfc_pgRNA_single_targeting_GI <- d.lfc_pgRNA_single_targeting %>%
  left_join(d.lfc_pgRNA_single_targeting_mean_lm_summary, by = "rep") %>%
  mutate(GI_score = single_target_CS - (intercept + slope * expected_CS)) 

# raw_GI = target_mean_CS - (control_intercept + control_slope * target_mean_expected_C
d.lfc_pgRNA_double_targeting_GI <- d.lfc_pgRNA_double_targeting %>%
  left_join(d.lfc_pgRNA_single_targeting_mean_lm_summary, by = "rep") %>%
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
## use lapply to run for each rep, then bind results
reps <- d.GI_scores_pgRNA_double %>%
  ungroup() %>%
  distinct(rep) %>% 
  pull()

results <- lapply(reps, function(i){
  ## get a vector of GI scores for all single-targeting ("control") pgRNAs for each rep
  single_GI_scores <- d.GI_scores_pgRNA_single %>%
    filter(rep == i) %>%
    pull(GI_score)
  
  ## get double-targeting pgRNAs for this rep, do a t-test to compare the double-
  ## targeting GI scores for each paralog pair to the control vector
  d.double_GI_scores <- d.GI_scores_pgRNA_double %>%
    filter(rep == i) %>%
    group_by(paralog_pair) %>%
    mutate(p_val = t.test(x = single_GI_scores,
                          y = GI_score,
                          paired = FALSE)$p.value) 
  
  ## adjust for multiple testing using the Benjamini-Hochberg method
  d.p_val <- d.double_GI_scores %>%
    dplyr::select(paralog_pair, p_val) %>%
    arrange(p_val) %>%
    distinct(p_val, .keep_all = TRUE) 
  
  p_vals <- d.p_val %>% 
    pull(p_val)
  
  fdr_vals <- p.adjust(p_vals, method = "BH")
  
  d.fdr <- tibble("fdr" = fdr_vals) %>%
    bind_cols(d.p_val) %>%
    dplyr::select(-p_val)
  
  ## add FDR values back into the double-targeting DF
  d.double_GI_scores <- left_join(d.double_GI_scores, d.fdr, by = "paralog_pair")
  
  return(d.double_GI_scores)
  
})

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
  group_by(rep, pgRNA_target) %>%
  mutate(mean_observed_CS = mean(observed_CS),
         mean_expected_CS = mean(expected_CS),
         mean_GI_score = mean(GI_score)) %>%
  distinct(pgRNA_target, .keep_all = TRUE) %>%
  dplyr::select(-c(pgRNA_id, contains("seq"), observed_CS, expected_CS, GI_score)) %>%
  ungroup()

d.GI_scores_target %>%
  group_by(rep) %>%
  summarize(n = n())

d.GI_scores_target_plot <- ggplot(d.GI_scores_target, aes(x = mean_expected_CS, y = mean_observed_CS, color = broad_target_type)) +
  geom_point() + 
  geom_abline(data = d.lfc_pgRNA_single_targeting_mean_lm_summary,
              aes(slope = slope, intercept = intercept),
              color = "gray55", size = 1) +
  labs(x = "expected_CRISPR_score", y = "observed_CRISPR_score") +
  plot_options +
  plot_theme +
  theme(aspect.ratio = square_ar) +
  facet_wrap(~rep)
save_plot(d.GI_scores_target_plot)

##Collapse target-level GI scores across reps - Added by Alice

d.GI_scores_target_rep <- d.GI_scores_target %>%
  group_by(pgRNA_target) %>%
  summarize(
    rep_collapse_mean_GI_score = mean(mean_GI_score),
    rep_collapse_mean_p_val = mean(p_val),
    rep_collapse_mean_fdr = mean(fdr),
    rep_collapse_mean_obs = mean(mean_observed_CS),
    rep_collapse_mean_exp = mean(mean_expected_CS)
    ) 

d.GI_scores_target_rep_collapsed <- d.GI_scores_target_rep %>%
  left_join(
    d.GI_scores_target %>% 
      select(paralog_pair, target_type, gene1_symbol, gene2_symbol, broad_target_type, pgRNA_target) %>%
      distinct(pgRNA_target, .keep_all = TRUE),   # one row per target
    by = "pgRNA_target"
  )

#####ADD MORE FIGURES HERE#######

## Save output
save_tbl(d.GI_scores_target)
save_tbl(d.GI_scores_pgRNA)




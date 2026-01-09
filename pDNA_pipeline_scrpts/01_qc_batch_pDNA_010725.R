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
library(pheatmap) # for making correlation heatmap
library(RColorBrewer) # for heatmap colors
library(kableExtra) # for formatting kables

args <- commandArgs(trailingOnly = TRUE)

base_dir  <- args[1]
cell_line <- args[2]



## get functions and vars from shared R script
source(file.path( base_dir, "00-shared_functions_and_variables.R"))

contour_palette <- colorRampPalette(brewer.pal(n = 9, name ="Spectral"))(50)


#Alice added lines to force directory creation

in_dir <- file.path( base_dir, "results", "pgRNA_counts")
if (!dir.exists(in_dir)) dir.create(in_dir, recursive = TRUE)

out_dir <- file.path( base_dir, "results", "pgRNA_counts_QC")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

d.counts <- read_tsv(file.path(in_dir, paste0(cell_line, "_counts_010526.txt")), col_names=TRUE)

## what does d.counts look like?
print_kbl(head(d.counts, 10))

## pivot_longer to tidy format
d.counts_long <- d.counts %>%
  dplyr::select(-c(seq_1, seq_2 )) %>%
  pivot_longer(!id, 
               names_to = "sample", 
               values_to = "count")

##Note 3800 is the number of pgRNAs in pgMI
##Alice comment - this step is to estimate average coverage of the library

d.coverage <- d.counts_long %>%
  group_by(sample) %>%
  summarize(sum = sum(count), 
            coverage = round((sum / 3800), 2))
print_kbl(d.coverage)
save_tbl(d.coverage)

coverage_histogram <- ggplot(d.coverage, aes(x = sample,
                                             y = coverage)) +
  geom_col() +
  geom_text(aes(label = round(coverage)),
            vjust = -0.5,
            size = 3) +
  theme_minimal() +
  labs(title = "Sample readcounts",
       x = "Sample",
       y = "Coverage (x)")

save_plot(coverage_histogram)

## convert sample to a factor 
d.counts_long <- d.counts_long %>%
  mutate(sample = factor(sample)) 

d.counts_cdf <- d.counts_long %>%
  group_by(sample) %>%
  mutate(count_norm = -log10((count+1)/sum(count)))

pgRNA_counts_cd <- ggplot(d.counts_cdf, aes(x = count_norm, color = sample)) +
  stat_ecdf() +
  labs(x = "-log10(count/total_count)", # bquote(~-log[10]~"(count/total_count)")
       y = "Expected_pgRNAs",
       color = "Sample") +  
  plot_options +
  plot_theme +
  theme(aspect.ratio = wide_ar)
pgRNA_counts_cd
save_plot(pgRNA_counts_cd)

## plot counts per million for each sample
d.counts_cpm <- d.counts_long %>%
  group_by(sample) %>%
  mutate(cpm = (((count)/sum(count))*1e6)) %>%
  mutate(log2_cpm = log2(cpm +1))

n_samples <- d.counts_long %>%
  distinct(sample) %>%
  nrow()

sample_cpm_histogram <- ggplot(d.counts_cpm, aes(x = log2_cpm, fill = sample)) +
  geom_histogram(color = "black", binwidth = 0.5) +
  plot_options +
  plot_theme +
  theme(aspect.ratio = wide_ar,
        legend.position = "none") +
  facet_wrap(~sample, scales = "free_y", ncol = ceiling(n_samples/2))
sample_cpm_histogram
save_plot(sample_cpm_histogram)

## sample correlation heatmap
d.counts_cpm_cor <- d.counts_cpm %>%
  dplyr::select(id, sample, cpm) %>%
  pivot_wider(names_from = "sample",
              values_from = "cpm") %>%
  dplyr::select(-id) %>%
  cor() %>%
  round(2) %>%
  data.frame()
print_kbl(d.counts_cpm_cor)

colors <- colorRampPalette(brewer.pal(n = 9, name ="YlGnBu"))(50)

sample_cor_heatmap_unfiltered <- pheatmap(d.counts_cpm_cor,
                                          col = colors,
                                          border_color = "white",
                                          cellwidth = 20, cellheight = 20,
                                          treeheight_row = 20, treeheight_col = 20,
                                          ## extra stuff
                                          cluster_rows = TRUE,
                                          cluster_cols = TRUE,
                                          cex = 1, clustering_distince_rows = "euclidean",
                                          cex = 1, clustering_distance_cols = "euclidean",
                                          cluster_method = "complete")
sample_cor_heatmap_unfiltered
save_plot(sample_cor_heatmap_unfiltered)

## flag pgRNAs with count = 0 at any time point
d.counts_cpm_filter <- d.counts_cpm %>%
  group_by(id) %>%
  mutate(zero_count = case_when(
    any(count == 0) ~ TRUE, ## if any value in the group = 0, set value to TRUE
    TRUE ~ FALSE)) %>% ## if above condition is not met, set value to FALSE
  ungroup()

## how many guides will be removed using this filter? (zero_count == TRUE)
d.summ <- d.counts_cpm_filter %>%
  dplyr::select(id, zero_count) %>%
  distinct(id, .keep_all = TRUE) %>%
  group_by(zero_count) %>%
  summarize(n = n()) %>%
  mutate(percent = round(((n/sum(n))*100), 2))
print_kbl(d.summ)

##Alice updated to use T0 counts
## flag pgRNAs with low T0 read counts

## filter pgRNA df for T0 reads only
d.counts_cpm_pDNA <- d.counts_cpm %>%
  filter(sample == "pDNA")

## what does the T0 read count distribution look like? 
pDNA_cpm_histogram <-  ggplot(d.counts_cpm_pDNA, aes(x = log2_cpm)) +
  geom_histogram(binwidth = 0.2,
                 color = "black", 
                 fill = "gray60") +
  plot_options +
  plot_theme +
  theme(aspect.ratio = wide_ar)
pDNA_cpm_histogram
save_plot(pDNA_cpm_histogram)

## determine T0 cutoff
d.counts_cpm_pDNA_stats <- d.counts_cpm_pDNA %>%
  summarize(median = median(log2_cpm),
            Q1 = quantile(log2_cpm, probs = 0.25),
            Q3 = quantile(log2_cpm, probs = 0.75),
            lower_outlier = (Q1 - 1.5*(Q3 - Q1)))
print_kbl(d.counts_cpm_pDNA_stats)

## save your selected cutoff as a variable
pDNA_cpm_cutoff <- d.counts_cpm_pDNA_stats %>%
  pull(lower_outlier) %>%
  unlist() %>%
  unname()

## This is now using T0 
## add your cutoff line to the T0 plot
pDNA_cpm_histogram_cutoff <- pDNA_cpm_histogram + 
  geom_vline(xintercept = pDNA_cpm_cutoff, ## adjust based on selected cutoff
             linetype = "dashed") 
pDNA_cpm_histogram_cutoff
save_plot(pDNA_cpm_histogram_cutoff)

## add filter variable to pgRNA df based on selected cutoff
d.counts_cpm_pDNA_filter <- d.counts_cpm_pDNA %>%
  mutate(low_pDNA_cpm = case_when(
    log2_cpm < pDNA_cpm_cutoff ~ TRUE, ## if plasmid log2_cpm < cutoff, set to TRUE
    TRUE ~ FALSE)) %>% ## if above condition is not met, set to FALSE
  ungroup()

## how many guides will be removed based on this filter?
d.summ <- d.counts_cpm_pDNA_filter %>%
  dplyr::select(id, low_pDNA_cpm) %>%
  distinct(id, .keep_all = TRUE) %>%
  group_by(low_pDNA_cpm) %>%
  summarize(n = n()) %>%
  mutate(percent = round(((n/sum(n))*100), 2))
print_kbl(d.summ)

## add low_T0_cpm variable to filter DF
d.counts_cpm_filter <- d.counts_cpm_pDNA_filter %>%
  dplyr::select(id, low_pDNA_cpm) %>%
  right_join(d.counts_cpm_filter, by = "id") %>%
  dplyr::select(id, sample:zero_count, low_pDNA_cpm) %>% ## reorder cols
  ungroup()

## how many pgRNAs will be removed by both filters?
d.counts_cpm_filter <- d.counts_cpm_filter %>%
  dplyr::select(id, zero_count, low_pDNA_cpm) %>%
  distinct(id, .keep_all = TRUE) %>%
  group_by(id) %>%
  mutate(rm_pgRNA = case_when(
    any(zero_count == TRUE | low_pDNA_cpm == TRUE) ~ TRUE,
    TRUE ~ FALSE)) %>%
  ungroup()

## write a function to do this summary? and return n(TRUE)?
d.summ <- d.counts_cpm_filter %>%
  group_by(rm_pgRNA) %>%
  summarize(n = n()) %>%
  mutate(percent = round(((n/sum(n))*100), 2))
print_kbl(d.summ) 

## confirm that n removed pgRNAs = # pgRNAs with low plasmid count + # pgRNAs with 0 reads at any time point
d.summ_complete <- d.counts_cpm_filter %>%
  mutate(flag_group = case_when(
    zero_count == TRUE & low_pDNA_cpm == TRUE ~ "both",
    zero_count == TRUE & low_pDNA_cpm == FALSE ~ "zero_count_only",
    low_pDNA_cpm == TRUE & zero_count == FALSE ~ "low_plasmid_cpm_only",
    low_pDNA_cpm == FALSE & zero_count == FALSE ~ "neither",
    TRUE ~ "error"
  )) %>%
  group_by(flag_group) %>%
  summarize(n = n())
print_kbl(d.summ_complete) 

d.counts_cpm_flag_long <- left_join(d.counts_cpm, d.counts_cpm_filter, by = "id")
save_tbl(d.counts_cpm_flag_long)

d.counts_cpm_flag_wide <- d.counts_cpm_flag_long %>%
  pivot_wider(names_from = sample, 
              values_from = count:log2_cpm,
              names_glue = "{sample}_{.value}")
save_tbl(d.counts_cpm_flag_wide)

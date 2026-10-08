## Project: Gut Microbiome in Spontaneous HIV-1 Controllers and Immunological Non Responders: a metagenomics study in 1,559 people with HIV
## Author: Ángela del Castillo Izquierdo
## Last time edited: 2026-10-07
##
## 4. Extra analysis -----------------------------------------------------------
## Loads working_objects.RData from script 0.
##
## Ratios use percent abundance and a pseudocount of 1: (a + 1) / (b + 1).
##   Prevotella / Bacteroides          genus table
##   Bacteroidales / Clostridiales     order table
## Group plots use a Wilcoxon test. The adjusted test is lm_phenotype
## (inverse-normal transform, then a linear model adjusted for age, MSM status,
## BMI, ethnicity, smoking status, recruitment center and log10 sequencing
## read count).
## Phylum differential abundance is EC vs NP on non-zero CLR, with those same
## covariates.
##
## Packages:
##   dplyr      1.2.0
##   ggplot2    4.0.2
##   patchwork  1.3.2
##   ggpubr     0.6.3
##   glmnet   4.1.10
##   reshape2 1.4.5
##   R.utils  2.13.0

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
})

root <- path.expand("~/Documents/Projects/HIV Elite Controllers/Files_Yue_2025_10_06/code_public")
source(file.path(root, "R", "R_functions.R"))
base::load(file.path(root, "data", "working_objects.RData"))
res_dir <- file.path(root, "results")
dir.create(res_dir, showWarnings = FALSE)

covar <- c("log10_counts", "MSM_derived", "CENTER", "AGE",
           "BMI_BASELINE", "SMOKING_CURRENT", "ETHNICITY")

abundance_ratio <- function(rel_mat, numerator, denominator) {
  x <- as.data.frame(rel_mat, check.names = FALSE) * 100
  stopifnot(numerator %in% names(x), denominator %in% names(x))
  ratio <- (x[[numerator]] + 1) / (x[[denominator]] + 1)
  stats::setNames(ratio, rownames(x))
}

prep_cov <- function(pheno_df) {
  d <- pheno_df[, covar, drop = FALSE]
  d$AGE <- as.integer(d$AGE)
  d$CENTER <- as.factor(as.character(d$CENTER))
  d$ETHNICITY <- as.factor(as.character(d$ETHNICITY))
  d$MSM_derived <- as.factor(d$MSM_derived)
  d$SMOKING_CURRENT <- as.factor(d$SMOKING_CURRENT)
  d[complete.cases(d), , drop = FALSE]
}

# One outcome column against one grouping factor.
run_outcome_lm <- function(y, pheno_df, group_col, cov_df, outcome_name) {
  pheno_df[[group_col]] <- factor(pheno_df[[group_col]])
  y_df <- data.frame(outcome = as.numeric(y), stringsAsFactors = FALSE)
  rownames(y_df) <- names(y)
  names(y_df) <- outcome_name
  samples <- Reduce(intersect, list(names(y), rownames(pheno_df), rownames(cov_df)))
  res <- lm_phenotype(
    y_df[samples, , drop = FALSE],
    pheno_df[samples, group_col, drop = FALSE],
    cov_df[samples, , drop = FALSE],
    colnames(cov_df)
  )
  data.frame(outcome = outcome_name,
             phenotype = group_col,
             Beta = res$Beta, SE = res$SE, P_value = res$p, N = res$N,
             stringsAsFactors = FALSE)
}

plot_ratio <- function(sample_ids, ratio, group, comparisons, ylab, title) {
  d <- data.frame(SampleID = sample_ids,
                  Ratio = ratio[sample_ids],
                  Group = group,
                  stringsAsFactors = FALSE)
  d <- d[complete.cases(d), ]
  d$Group <- factor(d$Group)
  ggplot(d, aes(x = Group, y = log2(Ratio))) +
    geom_jitter(width = 0.2, alpha = 0.7, colour = "grey") +
    geom_boxplot(width = 0.5, fill = "lightblue", alpha = 0.5, outliers = FALSE) +
    ggpubr::stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                               label = "p.signif", hide.ns = FALSE) +
    labs(x = "Control status", y = ylab, title = title) +
    theme_bw() +
    theme(axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, colour = "black"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())
}

# --- 4.1 Prevotella / Bacteroides --------------------------------------------
pb_ratio <- abundance_ratio(sp_re_metaphlan_genus, "g__Prevotella", "g__Bacteroides")

p_pb_hic <- plot_ratio(
  pheno$SampleID, pb_ratio, pheno$EC_lost,
  list(c("EC_persistent", "non_EC"), c("EC_transient", "non_EC")),
  "Prevotella/Bacteroides ratio (log2)",
  "HIC-persistent and HIC-transient vs NP"
)
p_pb_ec <- plot_ratio(
  pheno_EC$SampleID, pb_ratio, pheno_EC$EC_yesno,
  list(c("EC", "non_EC")),
  "Prevotella/Bacteroides ratio (log2)",
  "EC vs NP"
)
p_pb <- p_pb_hic / p_pb_ec

pb_lm <- run_outcome_lm(pb_ratio, pheno_EC, "EC_yesno", prep_cov(pheno_EC),
                        "Prevotella_Bacteroides")

# --- 4.2 Bacteroidales / Clostridiales ---------------------------------------
bc_ratio <- abundance_ratio(sp_re_metaphlan_order, "Bacteroidales", "Clostridiales")

p_bc_ec <- plot_ratio(
  pheno_EC$SampleID, bc_ratio, pheno_EC$EC_yesno,
  list(c("EC", "non_EC")),
  "Bacteroidales/Clostridiales ratio (log2)",
  "EC vs NP"
)
p_bc_hic <- plot_ratio(
  pheno$SampleID, bc_ratio, pheno$EC_lost,
  list(c("EC_persistent", "non_EC"), c("EC_transient", "non_EC")),
  "Bacteroidales/Clostridiales ratio (log2)",
  "HIC-persistent and HIC-transient vs NP"
)
p_bc <- p_bc_ec / p_bc_hic

bc_lm <- run_outcome_lm(bc_ratio, pheno_EC, "EC_yesno", prep_cov(pheno_EC),
                        "Bacteroidales_Clostridiales")

# --- 4.3 Phylum, EC vs NP -----------------------------------------------------
cov_phylum <- prep_cov(pheno_EC)
pheno_ec <- pheno_EC
pheno_ec$EC_yesno <- factor(pheno_ec$EC_yesno)
samples_p <- Reduce(intersect, list(
  rownames(sp_re_metaphlan_clr_phylum_non_zero),
  rownames(pheno_ec),
  rownames(cov_phylum)
))
phylum_res <- lm_phenotype(
  sp_re_metaphlan_clr_phylum_non_zero[samples_p, , drop = FALSE],
  pheno_ec[samples_p, "EC_yesno", drop = FALSE],
  cov_phylum[samples_p, , drop = FALSE],
  colnames(cov_phylum)
)
phylum_res <- phylum_res[!is.na(phylum_res$SE) & !is.nan(phylum_res$SE), ]
phylum_lm <- data.frame(
  Taxa = phylum_res$taxa,
  Phenotype = "EC vs NP",
  Beta = phylum_res$Beta,
  SE = phylum_res$SE,
  P_value = phylum_res$p,
  P_adjusted_FDR = phylum_res$fdr.p,
  N = phylum_res$N,
  stringsAsFactors = FALSE
)
phylum_lm <- phylum_lm[order(phylum_lm$P_adjusted_FDR), ]

ratio_lm <- bind_rows(pb_lm, bc_lm)
write.csv(ratio_lm, file.path(res_dir, "04_ratios_lm.csv"), row.names = FALSE)
write.csv(phylum_lm, file.path(res_dir, "04_phylum_lm_EC.csv"), row.names = FALSE)

ratio_lm
phylum_lm

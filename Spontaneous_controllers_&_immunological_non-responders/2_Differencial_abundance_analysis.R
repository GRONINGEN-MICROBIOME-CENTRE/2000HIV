## Project: Gut Microbiome in Spontaneous HIV-1 Controllers and Immunological Non Responders: a metagenomics study in 1,559 people with HIV
## Author: Ángela del Castillo Izquierdo
## Last time edited: 2026-10-08
##
## 2. Differential abundance ---------------------------------------------------
## Loads working_objects.RData from script 0.
##
## Linear models (lm_phenotype) use non-zero CLR abundance.
## Logistic models (lrm_phenotypes) use presence/absence.
## Both inverse-normal transform the phenotype and the covariates, then adjust
## for age, MSM status, BMI, ethnicity, smoking status, recruitment center and
## log10 sequencing read count. FDR is Benjamini-Hochberg across taxa in that
## model.
## Group columns are factor() with alphabetical levels, so the beta refers to
## that coding (EC = 1, non_EC = 2; Non-Responder = 1, Responder = 2;
## EC_persistent / EC_transient = 1, non_EC = 2).
##
## The same four contrasts are run for pathways. Pathways were not
## prevalence-filtered. Only EC vs NP, at genus and at SGB, is also fit with
## latest CD4+ T-cell count added to the covariates.
## Violin plots show taxa with FDR < 0.1 and the unadjusted lm p-value.
## The scatter highlights taxa with FDR < 0.1 in every contrast (filled) and
## taxa with FDR < 0.1 in at least one contrast (open).
##
## Packages:
##   dplyr      1.2.0
##   tidyr      1.3.2
##   ggplot2    4.0.2
##   patchwork  1.3.2
##   scales     1.4.0
##   glmnet     4.1.10
##   reshape2   1.4.5
##   R.utils    2.13.0

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

root <- path.expand("~/Documents/Projects/HIV Elite Controllers/Files_Yue_2025_10_06/code_public")
source(file.path(root, "R", "R_functions.R"))
base::load(file.path(root, "data", "working_objects.RData"))
res_dir <- file.path(root, "results")
dir.create(res_dir, showWarnings = FALSE)

covar <- c("log10_counts", "MSM_derived", "CENTER", "AGE",
           "BMI_BASELINE", "SMOKING_CURRENT", "ETHNICITY")
FDR_CUT <- 0.1

prep_cov <- function(pheno_df, extra = NULL) {
  d <- pheno_df[, c(covar, extra), drop = FALSE]
  d$AGE <- as.integer(d$AGE)
  d$CENTER <- as.factor(as.character(d$CENTER))
  d$ETHNICITY <- as.factor(as.character(d$ETHNICITY))
  d$MSM_derived <- as.factor(d$MSM_derived)
  d$SMOKING_CURRENT <- as.factor(d$SMOKING_CURRENT)
  if ("CD4_LATEST" %in% names(d)) d$CD4_LATEST <- as.numeric(d$CD4_LATEST)
  d[complete.cases(d), , drop = FALSE]
}

contrasts <- list(
  EC_vs_NP = list(pheno = pheno_EC, group = "EC_yesno",
                  cols = c(non_EC = "darkslategray3", EC = "salmon"),
                  levels = c("non_EC", "EC")),
  INR_vs_IR = list(pheno = pheno_IR, group = "INR_YESNO",
                   cols = c(`Non-Responder` = "plum", Responder = "#6666CC"),
                   levels = c("Non-Responder", "Responder")),
  HIC_persistent_vs_NP = list(pheno = pheno_persistent, group = "EC_lost",
                              cols = c(non_EC = "darkslategray3", EC_persistent = "sandybrown"),
                              levels = c("non_EC", "EC_persistent")),
  HIC_transient_vs_NP = list(pheno = pheno_transient, group = "EC_lost",
                             cols = c(non_EC = "darkslategray3", EC_transient = "lightgoldenrod"),
                             levels = c("non_EC", "EC_transient"))
)

run_lm <- function(y_mat, pheno_df, group_col, cov_df, contrast, level) {
  pheno_df[[group_col]] <- factor(pheno_df[[group_col]])
  samples <- Reduce(intersect, list(rownames(y_mat), rownames(pheno_df), rownames(cov_df)))
  res <- lm_phenotype(
    y_mat[samples, , drop = FALSE],
    pheno_df[samples, group_col, drop = FALSE],
    cov_df[samples, , drop = FALSE],
    colnames(cov_df)
  )
  res <- res[!is.na(res$SE) & !is.nan(res$SE), ]
  res$contrast <- contrast
  res$level <- level
  res$method <- "lm"
  res
}

run_lrm <- function(y_mat, pheno_df, group_col, cov_df, contrast, level) {
  pheno_df[[group_col]] <- factor(pheno_df[[group_col]])
  samples <- Reduce(intersect, list(rownames(y_mat), rownames(pheno_df), rownames(cov_df)))
  res <- lrm_phenotypes(
    y_mat[samples, , drop = FALSE],
    pheno_df[samples, group_col, drop = FALSE],
    cov_df[samples, , drop = FALSE],
    colnames(cov_df)
  )
  res <- res[!is.na(res$SE) & !is.nan(res$SE), ]
  res <- res[res$y_0_N > 0 & res$y_1_N > 0, ]
  res$taxa <- res$Y
  res$contrast <- contrast
  res$level <- level
  res$method <- "lrm"
  res
}

tidy_da <- function(res) {
  data.frame(Taxa = res$taxa,
             contrast = res$contrast,
             level = res$level,
             method = res$method,
             Beta = res$Beta,
             SE = res$SE,
             P_value = res$p,
             P_adjusted_FDR = res$fdr.p,
             N = res$N,
             stringsAsFactors = FALSE)
}

lm_mats <- list(SGB = sp_re_metaphlan_clr_non_zero,
                Genus = sp_re_metaphlan_clr_genus_non_zero,
                Pathway = re_path_clr_non_zero)
lrm_mats <- list(SGB = sp_re_metaphlan_binary,
                 Genus = sp_re_metaphlan_genus_binary,
                 Pathway = re_path_binary)
rel_mats <- list(SGB = sp_re_metaphlan, Genus = sp_re_metaphlan_genus)

lm_res <- list()
lrm_res <- list()
for (lv in names(lm_mats)) {
  for (nm in names(contrasts)) {
    sp <- contrasts[[nm]]
    cov_df <- prep_cov(sp$pheno)
    lm_res[[paste(lv, nm)]] <- run_lm(lm_mats[[lv]], sp$pheno, sp$group, cov_df, nm, lv)
    lrm_res[[paste(lv, nm)]] <- run_lrm(lrm_mats[[lv]], sp$pheno, sp$group, cov_df, nm, lv)
  }
}

# EC vs NP only, genus and SGB: same linear model with latest CD4+ T-cell count.
cd4_res <- list()
for (lv in c("SGB", "Genus")) {
  cov_cd4 <- prep_cov(pheno_EC, extra = "CD4_LATEST")
  cd4_res[[lv]] <- run_lm(lm_mats[[lv]], pheno_EC, "EC_yesno", cov_cd4,
                          "EC_vs_NP_plus_CD4", lv)
}

da_lm <- bind_rows(lapply(lm_res, tidy_da))
da_lrm <- bind_rows(lapply(lrm_res, tidy_da))
da_cd4 <- bind_rows(lapply(cd4_res, tidy_da))
write.csv(da_lm, file.path(res_dir, "02_DA_lm.csv"), row.names = FALSE)
write.csv(da_lrm, file.path(res_dir, "02_DA_lrm.csv"), row.names = FALSE)
write.csv(da_cd4, file.path(res_dir, "02_DA_lm_EC_plus_CD4.csv"), row.names = FALSE)
base::save(lm_res, lrm_res, cd4_res, da_lm, da_lrm, da_cd4,
     file = file.path(res_dir, "02_differential_abundance.RData"))

# --- Violin / box / jitter for lm hits, with the lm p-value ------------------
plot_hits <- function(res, clr_mat, pheno_df, group_col, cols, x_levels, title) {
  hits <- res$taxa[!is.na(res$fdr.p) & res$fdr.p < FDR_CUT]
  if (!length(hits)) return(NULL)
  samples <- intersect(rownames(clr_mat), rownames(pheno_df))
  long <- as.data.frame(clr_mat[samples, hits, drop = FALSE])
  long$SampleID <- rownames(long)
  long$group <- as.character(pheno_df[samples, group_col])
  long <- long %>%
    pivot_longer(all_of(hits), names_to = "taxa", values_to = "clr") %>%
    filter(is.finite(clr))
  long$group <- factor(long$group, levels = x_levels)
  ann <- res[res$taxa %in% hits, ]
  ann$label <- paste0("p = ", format.pval(ann$p, digits = 2, eps = 1e-3))
  ann$taxa <- factor(ann$taxa, levels = hits)
  long$taxa <- factor(long$taxa, levels = hits)
  ggplot(long, aes(x = group, y = clr, colour = group, fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.35, colour = "white", width = 0.8) +
    geom_jitter(width = 0.12, size = 1.1, alpha = 0.75) +
    geom_boxplot(outlier.shape = NA, width = 0.28, colour = "black", alpha = 0.4) +
    geom_text(data = ann, aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = 1.05, vjust = 1.4, size = 3) +
    facet_wrap(~ taxa, scales = "free_y") +
    scale_colour_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(x = NULL, y = "CLR abundance (non-zero)", title = title) +
    theme_bw() +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"))
}

hit_plots <- list()
for (key in names(lm_res)) {
  lv <- sub(" .*", "", key)
  nm <- sub("^\\S+ ", "", key)
  sp <- contrasts[[nm]]
  hit_plots[[key]] <- plot_hits(lm_res[[key]], lm_mats[[lv]], sp$pheno, sp$group,
                                sp$cols, sp$levels, paste(nm, lv, "lm FDR < 0.1"))
}
hit_plots <- hit_plots[!vapply(hit_plots, is.null, logical(1))]
p_lm_violins <- if (length(hit_plots)) wrap_plots(hit_plots, ncol = 1) else NULL

# Forest plot: EC vs NP taxa that are FDR < 0.1 without CD4, with and without CD4.
forest_one <- function(primary, with_cd4, level) {
  hits <- primary$taxa[!is.na(primary$fdr.p) & primary$fdr.p < FDR_CUT]
  if (!length(hits)) return(NULL)
  d <- bind_rows(
    primary %>% filter(taxa %in% hits) %>% mutate(model = "Without CD4"),
    with_cd4 %>% filter(taxa %in% hits) %>% mutate(model = "With CD4")
  ) %>%
    mutate(CI_lower = Beta - 1.96 * SE,
           CI_upper = Beta + 1.96 * SE,
           taxa = factor(taxa, levels = rev(hits)))
  ggplot(d, aes(x = Beta, y = taxa, colour = model)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                   height = 0.2, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5), size = 2.4) +
    scale_colour_manual(values = c("Without CD4" = "grey30", "With CD4" = "salmon")) +
    labs(x = "Effect size (beta)", y = NULL, colour = NULL, title = paste("EC vs NP,", level)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
}
forests <- lapply(names(cd4_res), function(lv) {
  key <- paste(lv, "EC_vs_NP")
  forest_one(lm_res[[key]], cd4_res[[lv]], lv)
})
forests <- forests[!vapply(forests, is.null, logical(1))]
p_forest <- if (length(forests)) wrap_plots(forests, ncol = 1) else NULL

# --- Stacked presence/absence bars for logistic-model hits -------------------
stack_one <- function(res, bin_mat, pheno_df, group_col, title) {
  hits <- unique(res$taxa[!is.na(res$fdr.p) & res$fdr.p < FDR_CUT])
  if (!length(hits)) return(NULL)
  samples <- intersect(rownames(bin_mat), rownames(pheno_df))
  d <- as.data.frame(bin_mat[samples, hits, drop = FALSE])
  d$group <- factor(pheno_df[samples, group_col])
  d <- d %>%
    pivot_longer(all_of(hits), names_to = "taxa", values_to = "present") %>%
    mutate(present = ifelse(present > 0, "Present", "Absent"))
  ggplot(d, aes(x = group, fill = present)) +
    geom_bar(position = "fill", width = 0.6, colour = "black") +
    facet_wrap(~ taxa, scales = "free_x") +
    scale_fill_manual(values = c(Absent = "white", Present = "black")) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Proportion of participants", fill = NULL, title = title) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"))
}
stack_plots <- list()
for (key in names(lrm_res)) {
  lv <- sub(" .*", "", key)
  nm <- sub("^\\S+ ", "", key)
  sp <- contrasts[[nm]]
  stack_plots[[key]] <- stack_one(lrm_res[[key]], lrm_mats[[lv]], sp$pheno, sp$group,
                                  paste(nm, lv, "lrm FDR < 0.1"))
}
stack_plots <- stack_plots[!vapply(stack_plots, is.null, logical(1))]
p_lrm_bars <- if (length(stack_plots)) wrap_plots(stack_plots, ncol = 1) else NULL

# --- Mean relative abundance x prevalence ------------------------------------
# Filled points: FDR < 0.1 in every contrast, in both the linear and the
# logistic model. Open red points: FDR < 0.1 in at least one of those models.
prevalence_scatter <- function(rel_mat, da_level, level_name) {
  x <- t(as.matrix(rel_mat))
  plot_data <- data.frame(
    Taxon = rownames(x),
    Mean_Rel_Abund = rowMeans(x, na.rm = TRUE) * 100,
    Prevalence = rowMeans(x > 0, na.rm = TRUE) * 100,
    stringsAsFactors = FALSE
  )
  sub <- da_level[da_level$level == level_name, ]
  n_contrasts <- length(unique(sub$contrast))
  hit_counts <- sub %>%
    filter(!is.na(P_adjusted_FDR), P_adjusted_FDR < FDR_CUT) %>%
    count(Taxa, method, name = "n_hit")
  any_hit <- unique(hit_counts$Taxa)
  both_methods <- hit_counts %>%
    group_by(Taxa) %>%
    summarise(n_method = n_distinct(method), n_hit = sum(n_hit), .groups = "drop") %>%
    filter(n_method == 2, n_hit >= 2 * n_contrasts)
  plot_data$Highlight <- ifelse(plot_data$Taxon %in% both_methods$Taxa, "Every contrast",
                         ifelse(plot_data$Taxon %in% any_hit, "Any contrast", "Other"))
  plot_data$Highlight <- factor(plot_data$Highlight,
                                levels = c("Other", "Any contrast", "Every contrast"))
  ggplot(plot_data, aes(x = Prevalence, y = Mean_Rel_Abund,
                        colour = Highlight, fill = Highlight, shape = Highlight)) +
    geom_point(alpha = 0.75, size = 2.2) +
    scale_y_log10() +
    scale_colour_manual(values = c(Other = "grey70", `Any contrast` = "red", `Every contrast` = "red")) +
    scale_fill_manual(values = c(Other = "grey70", `Any contrast` = NA, `Every contrast` = "red")) +
    scale_shape_manual(values = c(Other = 16, `Any contrast` = 21, `Every contrast` = 21)) +
    labs(x = "Prevalence (% of participants)",
         y = "Mean relative abundance (%)",
         title = level_name, colour = NULL, fill = NULL, shape = NULL) +
    theme_minimal() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.minor = element_blank())
}

da_all <- bind_rows(da_lm, da_lrm)
p_prev <- prevalence_scatter(sp_re_metaphlan_genus, da_all, "Genus") +
  prevalence_scatter(sp_re_metaphlan, da_all, "SGB")
p_prev

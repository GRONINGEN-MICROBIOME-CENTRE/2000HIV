## Project: Gut Microbiome in Spontaneous HIV-1 Controllers and Immunological Non Responders: a metagenomics study in 1,559 people with HIV
## Author: Ángela del Castillo Izquierdo
## Last time edited: 2026-10-07
##
## 1. Microbial profiling and diversity ---------------------------------------
## Loads working_objects.RData from script 0.
##
##   Characterised vs uncharacterised SGBs and genera, before and after the
##   >= 15 detection filter.
##   PERMANOVA of the main covariates, then of EC vs NP, INR vs IR,
##   HIC-persistent vs NP and HIC-transient vs NP. Distance is Euclidean on
##   the CLR table (taxa in columns). adonis2, by = "margin", 1,000 permutations.
##   Alpha diversity (observed richness, Shannon, Simpson 1-D) on relative
##   abundance. The test is lm_phenotype: inverse-normal transform, then a
##   linear model adjusted for age, MSM status, BMI, ethnicity, smoking status,
##   recruitment center and log10 sequencing read count.
##   Genus-level PCA for EC vs NP and INR vs IR, with 95% ellipses, segments
##   to the smaller group's centroid, and outlined centroids.
##
## Packages:
##   dplyr      1.2.0
##   tidyr      1.3.2
##   ggplot2    4.0.2
##   vegan      2.7.3
##   patchwork  1.3.2
##   R.utils    2.13.0
##   reshape2   1.4.5
##   glmnet     4.1.10

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(vegan)
  library(patchwork)
})

root <- path.expand("~/GitHub")
source(file.path(root, "R", "R_functions.R"))
base::load(file.path(root, "data", "working_objects.RData"))
res_dir <- file.path(root, "results")
dir.create(res_dir, showWarnings = FALSE)

covar <- c("log10_counts", "MSM_derived", "CENTER", "AGE",
           "BMI_BASELINE", "SMOKING_CURRENT", "ETHNICITY")
theme_div <- theme_bw() +
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", size = 14),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", colour = "black", size = 11))

# --- 1.1 Characterised / uncharacterised, before and after filtering ----------
# Uncharacterised SGB: species name starts with GGB.
# Uncharacterised genus: genus id starts with g__GGB.
status_counts <- function(ids, is_characterised, when) {
  data.frame(
    when = when,
    status = ifelse(is_characterised, "Characterised", "Uncharacterised"),
    stringsAsFactors = FALSE
  ) %>%
    count(when, status, name = "n")
}

sgb_before <- taxonomy_table$ID
sgb_char <- setNames(!startsWith(taxonomy_table$Species, "GGB"), taxonomy_table$ID)
genus_before <- genus_ids_unfiltered
genus_char <- setNames(!startsWith(genus_before, "g__GGB"), genus_before)

bar_df <- bind_rows(
  status_counts(sgb_before, sgb_char[sgb_before], "Before") %>% mutate(level = "SGB"),
  status_counts(colnames(sp_re_metaphlan), sgb_char[colnames(sp_re_metaphlan)], "After") %>%
    mutate(level = "SGB"),
  status_counts(genus_before, genus_char[genus_before], "Before") %>% mutate(level = "Genus"),
  status_counts(colnames(sp_re_metaphlan_genus), genus_char[colnames(sp_re_metaphlan_genus)], "After") %>%
    mutate(level = "Genus")
)
bar_df$when <- factor(bar_df$when, levels = c("Before", "After"))
bar_df$status <- factor(bar_df$status, levels = c("Characterised", "Uncharacterised"))

p_char <- ggplot(bar_df, aes(x = when, y = n, fill = status)) +
  geom_col(width = 0.7, colour = "black") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.5) +
  facet_wrap(~ level) +
  scale_fill_manual(values = c("Characterised" = "#1982c4", "Uncharacterised" = "#ffca3a")) +
  labs(x = NULL, y = "Number of taxa", fill = NULL,
       title = "Characterised and uncharacterised taxa before and after the prevalence filter") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
# --- 1.2 PERMANOVA ------------------------------------------------------------
run_permanova <- function(clr_mat, pheno_df, group_col = NULL) {
  use <- intersect(rownames(clr_mat), rownames(pheno_df))
  meta <- pheno_df[use, c(covar, group_col), drop = FALSE]
  meta <- meta[complete.cases(meta), , drop = FALSE]
  for (v in intersect(c("CENTER", "SMOKING_CURRENT", "MSM_derived", "ETHNICITY", group_col), names(meta))) {
    meta[[v]] <- as.factor(meta[[v]])
  }
  x <- clr_mat[rownames(meta), , drop = FALSE]
  dist_mat <- as.matrix(vegan::vegdist(x, method = "euclidean"))
  rhs <- paste(c(group_col, covar), collapse = " + ")
  set.seed(1)
  fit <- vegan::adonis2(as.formula(paste("dist_mat ~", rhs)),
                        data = meta[rownames(dist_mat), , drop = FALSE],
                        permutations = 1000, by = "margin")
  out <- as.data.frame(fit)
  out$Variable <- rownames(out)
  out$R2_percent <- out$R2 * 100
  out
}

permanova_bar <- function(fit_df, title) {
  d <- fit_df[!fit_df$Variable %in% c("Residual", "Total"), ]
  d <- d[order(d$R2_percent, decreasing = TRUE), ]
  d$Variable <- factor(d$Variable, levels = d$Variable)
  ggplot(d, aes(x = Variable, y = R2_percent)) +
    geom_col(fill = "skyblue", colour = "black", width = 0.7) +
    labs(x = NULL, y = "Beta-diversity variance explained (%)", title = title) +
    theme_minimal() +
    theme(axis.text.x = element_text(colour = "black", angle = 30, hjust = 1),
          axis.text.y = element_text(colour = "black"),
          axis.line = element_line(colour = "black"),
          panel.grid.minor = element_blank())
}

levels_fit <- list(SGB = sp_re_metaphlan_clr, Genus = sp_re_metaphlan_clr_genus)
covariate_fits <- lapply(levels_fit, run_permanova, pheno_df = pheno)
p_cov <- permanova_bar(covariate_fits$SGB, "SGB, covariates") +
  permanova_bar(covariate_fits$Genus, "Genus, covariates")
group_specs <- list(
  EC_vs_NP = list(pheno = pheno_EC, group = "EC_yesno"),
  INR_vs_IR = list(pheno = pheno_IR, group = "INR_YESNO"),
  HIC_persistent_vs_NP = list(pheno = pheno_persistent, group = "EC_lost"),
  HIC_transient_vs_NP = list(pheno = pheno_transient, group = "EC_lost")
)
group_fits <- lapply(names(levels_fit), function(lv) {
  lapply(names(group_specs), function(nm) {
    sp <- group_specs[[nm]]
    fit <- run_permanova(levels_fit[[lv]], sp$pheno, sp$group)
    fit$level <- lv
    fit$contrast <- nm
    fit
  })
})
names(group_fits) <- names(levels_fit)
group_fit_df <- bind_rows(lapply(group_fits, bind_rows))
write.csv(group_fit_df, file.path(res_dir, "01_permanova_groups.csv"), row.names = FALSE)
write.csv(bind_rows(covariate_fits, .id = "level"),
          file.path(res_dir, "01_permanova_covariates.csv"), row.names = FALSE)

# R2 of the grouping term only, for a compact bar.
group_r2 <- group_fit_df %>%
  filter(Variable %in% c("EC_yesno", "INR_YESNO", "EC_lost")) %>%
  mutate(contrast = factor(contrast, levels = names(group_specs)))
p_group <- ggplot(group_r2, aes(x = contrast, y = R2_percent, fill = level)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, colour = "black") +
  labs(x = NULL, y = "Beta-diversity variance explained (%)", fill = NULL,
       title = "PERMANOVA, grouping term (marginal)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"))
# --- 1.3 Alpha diversity ------------------------------------------------------
alpha_table <- function(rel_mat, pheno_df) {
  use <- intersect(rownames(rel_mat), rownames(pheno_df))
  x <- rel_mat[use, , drop = FALSE]
  data.frame(SampleID = use,
             Shannon = vegan::diversity(x, index = "shannon"),
             Simpson = vegan::diversity(x, index = "simpson"),
             Observed_Richness = vegan::specnumber(x),
             stringsAsFactors = FALSE)
}

# factor() uses alphabetical levels. as.numeric() inside lm_phenotype codes the
# first level as 1. Beta is the slope for that coding after the inverse-normal
# transform.
run_alpha_lm <- function(alpha_df, pheno_df, group_col, label, level) {
  merged <- merge(alpha_df, pheno_df, by = "SampleID")
  rownames(merged) <- merged$SampleID
  merged[[group_col]] <- factor(merged[[group_col]])
  merged$AGE <- as.integer(merged$AGE)
  merged$CENTER <- as.factor(as.character(merged$CENTER))
  merged$ETHNICITY <- as.factor(as.character(merged$ETHNICITY))
  merged$MSM_derived <- as.factor(merged$MSM_derived)
  merged$SMOKING_CURRENT <- as.factor(merged$SMOKING_CURRENT)
  cov_df <- merged[, covar, drop = FALSE]
  cov_df <- cov_df[complete.cases(cov_df), , drop = FALSE]
  samples <- intersect(rownames(merged), rownames(cov_df))
  y <- merged[samples, c("Shannon", "Simpson", "Observed_Richness"), drop = FALSE]
  res <- lm_phenotype(y, merged[samples, group_col, drop = FALSE], cov_df[samples, , drop = FALSE], covar)
  data.frame(contrast = label,
             level = level,
             metric = res$taxa,
             Beta = res$Beta,
             SE = res$SE,
             P_value = res$p,
             N = res$N,
             stringsAsFactors = FALSE)
}

plot_alpha <- function(alpha_df, pheno_df, group_col, cols, x_levels, lm_df, title) {
  merged <- merge(alpha_df, pheno_df[, c("SampleID", group_col)], by = "SampleID")
  long <- merged %>%
    pivot_longer(c(Shannon, Simpson, Observed_Richness), names_to = "Metric", values_to = "Value")
  long$group <- factor(long[[group_col]], levels = x_levels)
  ann <- lm_df %>%
    transmute(Metric = metric,
              label = paste0("lm p = ", format.pval(P_value, digits = 2, eps = 1e-3)))
  ggplot(long, aes(x = group, y = Value, colour = group, fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.3, colour = "white", width = 0.8) +
    geom_jitter(width = 0.13, size = 1.1, alpha = 0.8) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.35, colour = "black") +
    geom_text(data = ann, aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = 1.1, vjust = 1.4, size = 3.2) +
    facet_wrap(~ Metric, scales = "free_y") +
    scale_colour_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(y = "Alpha diversity", title = title) +
    theme_div +
    theme(legend.position = "none")
}

alpha_specs <- list(
  EC_vs_NP = list(pheno = pheno_EC, group = "EC_yesno",
                  levels = c("non_EC", "EC"),
                  cols = c(non_EC = "darkslategray3", EC = "salmon")),
  INR_vs_IR = list(pheno = pheno_IR, group = "INR_YESNO",
                   levels = c("Non-Responder", "Responder"),
                   cols = c(`Non-Responder` = "plum", Responder = "#6666CC")),
  HIC_persistent_vs_NP = list(pheno = pheno_persistent, group = "EC_lost",
                              levels = c("non_EC", "EC_persistent"),
                              cols = c(non_EC = "darkslategray3", EC_persistent = "sandybrown")),
  HIC_transient_vs_NP = list(pheno = pheno_transient, group = "EC_lost",
                             levels = c("non_EC", "EC_transient"),
                             cols = c(non_EC = "darkslategray3", EC_transient = "lightgoldenrod"))
)
rel_by_level <- list(SGB = sp_re_metaphlan, Genus = sp_re_metaphlan_genus)

alpha_results <- list()
alpha_plots <- list()
for (lv in names(rel_by_level)) {
  for (nm in names(alpha_specs)) {
    sp <- alpha_specs[[nm]]
    a <- alpha_table(rel_by_level[[lv]], sp$pheno)
    fit <- run_alpha_lm(a, sp$pheno, sp$group, nm, lv)
    alpha_results[[paste(lv, nm, sep = "__")]] <- fit
    alpha_plots[[paste(lv, nm, sep = "__")]] <- plot_alpha(
      a, sp$pheno, sp$group, sp$cols, sp$levels, fit,
      paste(nm, lv)
    )
  }
}
alpha_results_df <- bind_rows(alpha_results)
write.csv(alpha_results_df, file.path(res_dir, "01_alpha_diversity_lm.csv"), row.names = FALSE)

p_alpha <- wrap_plots(alpha_plots, ncol = 2)
# --- 1.4 Genus PCA: EC vs NP and INR vs IR ------------------------------------
pca_spider <- function(clr_mat, pheno_df, group_col, group_map, cols, spider_group) {
  use <- intersect(rownames(clr_mat), rownames(pheno_df))
  meta <- pheno_df[use, c(covar, group_col), drop = FALSE]
  meta <- meta[complete.cases(meta[, covar, drop = FALSE]), , drop = FALSE]
  x <- clr_mat[rownames(meta), , drop = FALSE]
  fit <- prcomp(x)
  pca_df <- as.data.frame(fit$x)
  raw <- as.character(meta[rownames(pca_df), group_col])
  pca_df$group <- unname(group_map[raw])
  pca_df <- pca_df[!is.na(pca_df$group), ]
  pca_df$group <- factor(pca_df$group, levels = names(cols))
  cents <- aggregate(cbind(PC1, PC2) ~ group, data = pca_df, FUN = mean)
  names(cents)[names(cents) %in% c("PC1", "PC2")] <- c("PC1_cent", "PC2_cent")
  pca_df <- merge(pca_df, cents, by = "group", all.x = TRUE, sort = FALSE)
  spider <- pca_df[pca_df$group == spider_group, ]
  var_pct <- round(100 * summary(fit)$importance[2, 1:2], 1)
  ggplot(pca_df, aes(x = PC1, y = PC2, colour = group, fill = group)) +
    stat_ellipse(geom = "polygon", linetype = 0, level = 0.95, alpha = 0.22, show.legend = FALSE) +
    geom_segment(data = spider,
                 aes(x = PC1, y = PC2, xend = PC1_cent, yend = PC2_cent, colour = group),
                 linewidth = 0.28, alpha = 0.8, show.legend = FALSE) +
    geom_point(size = 1.2, alpha = 0.9) +
    geom_point(data = cents, aes(x = PC1_cent, y = PC2_cent, fill = group),
               inherit.aes = FALSE, shape = 21, size = 2.8, colour = "black", stroke = 0.7,
               show.legend = FALSE) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_fill_manual(values = cols, name = NULL) +
    labs(x = paste0("PC1 (", var_pct[1], "% variance)"),
         y = paste0("PC2 (", var_pct[2], "% variance)"),
         title = paste0("N = ", format(nrow(pca_df), big.mark = ","))) +
    theme_bw(base_size = 12) +
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.background = element_blank(), panel.grid = element_blank())
}

p_pca_ec <- pca_spider(
  sp_re_metaphlan_clr_genus, pheno_EC, "EC_yesno",
  c(EC = "EC", non_EC = "NP"),
  c(EC = "#E07A5F", NP = "#7EBDC2"),
  "EC"
)
p_pca_inr <- pca_spider(
  sp_re_metaphlan_clr_genus, pheno_IR, "INR_YESNO",
  c(`Non-Responder` = "INR", Responder = "IR"),
  c(INR = "#E889A0", IR = "#5B5EA6"),
  "INR"
)
p_pca <- p_pca_ec + p_pca_inr + plot_annotation(tag_levels = "A")
p_char
p_cov
p_group
p_alpha
p_pca

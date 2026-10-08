## Project: Gut Microbiome in Spontaneous HIV-1 Controllers and Immunological Non Responders: a metagenomics study in 1,559 people with HIV
## Author: Ángela del Castillo Izquierdo
## Last time edited: 2026-10-07
##
## 3. Link taxa to phenotypes --------------------------------------------------
## Loads working_objects.RData from script 0.
## Absolute CD4 recovery is the latest CD4 count minus the CD4 nadir.
## Relative CD4 recovery is the latest CD4 count divided by the CD4 nadir.
##
## Panel A: EC vs NP and INR vs IR, with mirrored prevalence bars and delta CLR
##   above the heatmaps.
## Panel B: the same taxa in immune Responders only.
## One black star when FDR < 0.1, on the heatmap, the prevalence bars and the
##   delta-CLR markers. No star otherwise.
##
## Packages:
##   dplyr      1.2.0
##   tidyr      1.3.2
##   ggplot2    4.0.2
##   patchwork  1.3.2
##   openxlsx   4.2.8.1
##   scales     1.4.0
##   glmnet     4.1.10
##   reshape2   1.4.5
##   R.utils    2.13.0

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

root <- path.expand("~/Documents/Projects/HIV Elite Controllers/Files_Yue_2025_10_06/code_public")
source(file.path(root, "R", "R_functions.R"))
base::load(file.path(root, "data", "working_objects.RData"))
proj_dir <- file.path(root, "results")
dir.create(proj_dir, recursive = TRUE, showWarnings = FALSE)

## 1. Phenotypes ---------------------------------------------------------------
# Absolute recovery = CD4_LATEST - CD4_NADIR
# Relative recovery = CD4_LATEST / CD4_NADIR
add_recovery <- function(pheno_df) {
  nadir <- suppressWarnings(as.numeric(pheno_df$CD4_NADIR))
  latest <- suppressWarnings(as.numeric(pheno_df$CD4_LATEST))
  pheno_df$CD4_RECOVERY <- latest - nadir
  ratio <- latest / nadir
  ratio[!is.finite(ratio)] <- NA
  pheno_df$CD4_RECOVERY_RATIO <- ratio
  pheno_df
}

recode_cdc_stage <- function(pheno_df) {
  if (!"HIV_STAGE_CDC" %in% colnames(pheno_df)) return(pheno_df)
  # HIV_STAGE_CDC codebook: 0 = Stage 1, 1 = Stage 2, 2 = Stage 3, 3 = Unknown.
  # Recode only if still 0-based (skip if already shifted to 1/2/3).
  if (any(pheno_df$HIV_STAGE_CDC == 0, na.rm = TRUE)) {
    pheno_df$HIV_STAGE_CDC[pheno_df$HIV_STAGE_CDC == 3] <- NA
    pheno_df$HIV_STAGE_CDC <- pheno_df$HIV_STAGE_CDC + 1
  }
  pheno_df
}

pheno_EC <- recode_cdc_stage(add_recovery(as.data.frame(pheno_EC)))
pheno_IR <- recode_cdc_stage(add_recovery(as.data.frame(pheno_IR)))
rownames(pheno_EC) <- as.character(pheno_EC$SampleID)
rownames(pheno_IR) <- as.character(pheno_IR$SampleID)

# Viral load is pre-ART. CD4, CD8 and the CD4/CD8 ratio are the latest counts.
drug_vars <- c(
  "CD4_NADIR", "VL_PRECART",
  "CD4_LATEST", "CD8_LATEST", "CD4CD8_LATEST",
  "HIV_DURATION", "CART_DURATION", "HIV_STAGE_CDC",
  "CD4_RECOVERY", "CD4_RECOVERY_RATIO"
)
drug_vars <- intersect(drug_vars, colnames(pheno_IR))
pheno_labels <- c(
  CD4_NADIR = "CD4+ T cell nadir",
  VL_PRECART = "Viral load (copies/mL)",
  CD4_LATEST = "CD4+ T cell counts",
  CD8_LATEST = "CD8+ T cell counts",
  CD4CD8_LATEST = "CD4+/CD8+",
  HIV_DURATION = "HIV duration (years)",
  CART_DURATION = "ART duration (years)",
  HIV_STAGE_CDC = "CDC stage",
  CD4_RECOVERY = "Absolute CD4+ T cell recovery",
  CD4_RECOVERY_RATIO = "Relative CD4+ T cell recovery"
)

## 2. Covariates ---------------------------------------------------------------
prep_complete <- function(pheno_df) {
  pheno_df$AGE <- as.integer(pheno_df$AGE)
  pheno_df$CENTER <- as.factor(as.character(pheno_df$CENTER))
  pheno_df$ETHNICITY <- as.factor(as.character(pheno_df$ETHNICITY))
  pheno_df$MSM_derived <- as.factor(pheno_df$MSM_derived)
  pheno_df$SMOKING_CURRENT <- as.factor(pheno_df$SMOKING_CURRENT)
  d <- data.frame(
    log10_counts     = pheno_df$log10_counts,
    MSM_derived      = pheno_df$MSM_derived,
    ETHNICITY        = pheno_df$ETHNICITY,
    AGE              = pheno_df$AGE,
    SMOKING_CURRENT  = pheno_df$SMOKING_CURRENT,
    BMI_BASELINE     = pheno_df$BMI_BASELINE,
    CENTER           = pheno_df$CENTER
  )
  rownames(d) <- rownames(pheno_df)
  d[complete.cases(d), , drop = FALSE]
}

pheno_complete_ec <- prep_complete(pheno_EC)
pheno_complete_ir <- prep_complete(pheno_IR)
covar_complete <- colnames(pheno_complete_ir)

pheno_IR_IR <- pheno_IR[pheno_IR$INR_YESNO %in% "Responder", ]


## 3. Taxon sets ---------------------------------------------------------------
# Taxa with FDR < 0.1 in the earlier EC vs NP or INR vs IR comparisons.
taxa_ec_lm_genus  <- c("g__GGB3005", "g__GGB9172")
taxa_inr_lm_genus <- c("g__GGB269", "g__GGB9690", "g__GGB1618")
taxa_inr_lm_sgb   <- c("t__SGB376", "t__SGB15195", "t__SGB15318")
taxa_inr_lrm_sgb  <- c("t__SGB14754")

taxa_lm_genus <- c(taxa_ec_lm_genus, taxa_inr_lm_genus)
taxa_lm_sgb   <- taxa_inr_lm_sgb
taxa_inr_all  <- c(taxa_inr_lm_genus, taxa_inr_lm_sgb, taxa_inr_lrm_sgb)
taxa_order    <- c(taxa_ec_lm_genus, taxa_inr_all)

stopifnot(all(taxa_lm_genus %in% colnames(sp_re_metaphlan_clr_genus_non_zero)))
stopifnot(all(taxa_lm_sgb   %in% colnames(sp_re_metaphlan_clr_non_zero)))
stopifnot(all(taxa_inr_lrm_sgb %in% colnames(sp_re_metaphlan_binary)))

taxon_labs <- c(
  g__GGB3005  = "GGB3005\nFirmicutes",
  g__GGB9172  = "GGB9172\nFirmicutes",
  g__GGB269   = "GGB269\nMethanomassiliicoccaceae",
  g__GGB9690  = "GGB9690\nFirmicutes",
  g__GGB1618  = "GGB1618\nBacteroidales",
  t__SGB376   = "SGB376\nMethanomassiliicoccaceae",
  t__SGB15195 = "SGB15195\nFirmicutes",
  t__SGB15318 = "Faecalibacterium prausnitzii\nSGB15318",
  t__SGB14754 = "Bacteroides\nSGB14754"
)
axis_cols <- c(
  setNames(rep("salmon", length(taxa_ec_lm_genus)), taxa_ec_lm_genus),
  setNames(rep("plum",   length(taxa_inr_all)),     taxa_inr_all)
)

annotate_df <- function(combined, method, significant_in, level, cohort) {
  taxa_col <- if ("taxa" %in% names(combined)) combined$taxa else combined$Y
  pheno_col <- if ("phenotype" %in% names(combined)) combined$phenotype else combined$X
  data.frame(
    Taxa            = taxa_col,
    Phenotype       = pheno_col,
    Method          = method,
    Significant_in  = significant_in,
    Level           = level,
    Cohort          = cohort,
    Beta            = combined$Beta,
    SE              = combined$SE,
    P_value         = combined$p,
    P_adjusted_FDR  = combined$fdr.p,
    N               = combined$N,
    stringsAsFactors = FALSE
  )
}

run_lm_block <- function(y_mat, pheno_df, cov_df, samples, taxa, phenotypes) {
  result_list <- list()
  for (drug in phenotypes) {
    result_list[[drug]] <- lm_phenotype(
      y_mat[samples, taxa, drop = FALSE],
      pheno_df[samples, drug, drop = FALSE],
      cov_df[samples, , drop = FALSE],
      colnames(cov_df)
    )
  }
  bind_rows(result_list, .id = "model")
}

run_lrm_block <- function(y_mat, pheno_df, cov_df, samples, taxa, phenotypes) {
  result_list <- list()
  for (drug in phenotypes) {
    result_list[[drug]] <- lrm_phenotypes(
      y_mat[samples, taxa, drop = FALSE],
      pheno_df[samples, drug, drop = FALSE],
      cov_df[samples, , drop = FALSE],
      colnames(cov_df)
    )
  }
  out <- bind_rows(result_list, .id = "model")
  out <- out[out$Y %in% taxa, ]
  out$fdr.p <- p.adjust(out$p, method = "fdr")
  out$taxa <- out$Y
  out
}

## 4. Panel A: EC taxa in pheno_EC --------------------------------------------
sample_ec_genus <- Reduce(intersect, list(
  rownames(pheno_EC),
  rownames(pheno_complete_ec),
  rownames(sp_re_metaphlan_clr_genus_non_zero)
))
table(pheno_EC[sample_ec_genus, ]$EC_yesno)

combined_ec_gen <- run_lm_block(
  sp_re_metaphlan_clr_genus_non_zero, pheno_EC, pheno_complete_ec,
  sample_ec_genus, taxa_ec_lm_genus, drug_vars
)
format_df_ec <- annotate_df(combined_ec_gen, "LM", "EC-NP", "genus", "EC-NP")

## 5. Panel A: INR taxa in pheno_IR -------------------------------------------
sample_ir_genus <- Reduce(intersect, list(
  rownames(pheno_IR),
  rownames(pheno_complete_ir),
  rownames(sp_re_metaphlan_clr_genus_non_zero)
))
sample_ir_sgb <- Reduce(intersect, list(
  rownames(pheno_IR),
  rownames(pheno_complete_ir),
  rownames(sp_re_metaphlan_clr_non_zero)
))
sample_ir_lrm <- Reduce(intersect, list(
  rownames(pheno_IR),
  rownames(pheno_complete_ir),
  rownames(sp_re_metaphlan_binary)
))
table(pheno_IR[sample_ir_genus, ]$INR_YESNO)

combined_ir_gen <- run_lm_block(
  sp_re_metaphlan_clr_genus_non_zero, pheno_IR, pheno_complete_ir,
  sample_ir_genus, taxa_inr_lm_genus, drug_vars
)
combined_ir_sgb <- run_lm_block(
  sp_re_metaphlan_clr_non_zero, pheno_IR, pheno_complete_ir,
  sample_ir_sgb, taxa_inr_lm_sgb, drug_vars
)
combined_ir_lrm <- run_lrm_block(
  sp_re_metaphlan_binary, pheno_IR, pheno_complete_ir,
  sample_ir_lrm, taxa_inr_lrm_sgb, drug_vars
)

format_df_inr <- rbind(
  annotate_df(combined_ir_gen, "LM",  "INR-IR", "genus", "INR-IR"),
  annotate_df(combined_ir_sgb, "LM",  "INR-IR", "SGB",   "INR-IR"),
  annotate_df(combined_ir_lrm, "LRM", "INR-IR", "SGB",   "INR-IR")
)

## 6. Panel B: all taxa in Responders only ------------------------------------
sample_genus <- Reduce(intersect, list(
  pheno_IR_IR$SampleID,
  rownames(pheno_complete_ir),
  rownames(sp_re_metaphlan_clr_genus_non_zero)
))
sample_sgb <- Reduce(intersect, list(
  pheno_IR_IR$SampleID,
  rownames(pheno_complete_ir),
  rownames(sp_re_metaphlan_clr_non_zero)
))
sample_lrm <- Reduce(intersect, list(
  pheno_IR_IR$SampleID,
  rownames(pheno_complete_ir),
  rownames(sp_re_metaphlan_binary)
))
table(pheno_IR[sample_genus, ]$INR_YESNO)

combined_lm_genus <- run_lm_block(
  sp_re_metaphlan_clr_genus_non_zero, pheno_IR, pheno_complete_ir,
  sample_genus, taxa_lm_genus, drug_vars
)
combined_lm_sgb <- run_lm_block(
  sp_re_metaphlan_clr_non_zero, pheno_IR, pheno_complete_ir,
  sample_sgb, taxa_lm_sgb, drug_vars
)
combined_lrm <- run_lrm_block(
  sp_re_metaphlan_binary, pheno_IR, pheno_complete_ir,
  sample_lrm, taxa_inr_lrm_sgb, drug_vars
)

format_df_lm_genus <- annotate_df(combined_lm_genus, "LM",
  ifelse(combined_lm_genus$taxa %in% taxa_ec_lm_genus, "EC-NP", "INR-IR"),
  "genus", "IR")
format_df_lm_sgb <- annotate_df(combined_lm_sgb, "LM", "INR-IR", "SGB", "IR")
format_df_lrm    <- annotate_df(combined_lrm,    "LRM", "INR-IR", "SGB", "IR")

## 7. Combined tables ---------------------------------------------------------
format_df_within <- rbind(format_df_ec, format_df_inr)
format_df <- rbind(format_df_lm_genus, format_df_lm_sgb, format_df_lrm)

finish_table <- function(df) {
  df$Taxa <- factor(df$Taxa, levels = taxa_order)
  df$Phenotype <- factor(df$Phenotype, levels = drug_vars)
  df <- df[order(df$Phenotype, df$Taxa), ]
  rownames(df) <- NULL
  df
}
format_df_within <- finish_table(format_df_within)
format_df <- finish_table(format_df)

xlsx_out <- file.path(proj_dir, "03_linked_taxa_phenotypes.xlsx")
openxlsx::write.xlsx(
  list(`Within study groups` = format_df_within, Responders = format_df),
  file = xlsx_out
)

## 8. Two-panel heatmap (A within-cohort, B Responders) -----------------------
# One black star if FDR < 0.1, on the heatmap and on the prevalence and
# delta-CLR markers. No star otherwise.
annotate_sig <- function(p, fdr.p) {
  hit <- !is.na(fdr.p) & fdr.p < 0.1
  list(
    signif = ifelse(hit, "*", ""),
    sig_level = ifelse(hit, "FDR", "None")
  )
}

prep_heat <- function(df) {
  df %>%
    dplyr::mutate(
      signif = ifelse(!is.na(P_adjusted_FDR) & P_adjusted_FDR < 0.1, "*", ""),
      sig_level = ifelse(signif == "*", "FDR", "None"),
      taxa_lab = factor(unname(taxon_labs[as.character(Taxa)]),
                        levels = unname(taxon_labs[taxa_order]))
    )
}

heat_within <- prep_heat(format_df_within)
heat_ir     <- prep_heat(format_df)
beta_lim <- max(
  0.7,
  quantile(abs(c(heat_within$Beta, heat_ir$Beta)), 0.95, na.rm = TRUE),
  na.rm = TRUE
)

n_lab <- function(n) paste0("N = ", format(n, big.mark = ",", scientific = FALSE))

fill_scale <- scale_fill_gradient2(
  low = "blue", mid = "white", high = "red",
  midpoint = 0, limits = c(-beta_lim, beta_lim), oob = scales::squish
)
sig_scale <- scale_color_manual(
  values = c(FDR = "black", None = "transparent"),
  breaks = "FDR",
  name = "FDR < 0.1"
)

cols_ec <- c("EC" = "salmon", "NP" = "darkslategray3")
cols_ir <- c("INR" = "plum", "IR" = "#6666CC")

clr_mat_for <- function(tx) {
  if (startsWith(tx, "g__")) sp_re_metaphlan_clr_genus_non_zero else sp_re_metaphlan_clr_non_zero
}
bin_mat_for <- function(tx) {
  if (startsWith(tx, "g__")) sp_re_metaphlan_genus_binary else sp_re_metaphlan_binary
}

# Presence % per group; g1 plotted upward, g2 downward.
# Stars are attached later from the genome-wide LRM (presence ~ EC_yesno / INR_YESNO).
taxon_prev <- function(taxa, samples, group_vec, g1, g2) {
  out <- lapply(taxa, function(tx) {
    empty <- data.frame(
      taxa = tx, group = factor(c(g1, g2), levels = c(g1, g2)),
      prevalence = NA_real_, signed = NA_real_,
      stringsAsFactors = FALSE
    )
    bin <- bin_mat_for(tx)
    if (!tx %in% colnames(bin)) return(empty)
    s <- intersect(samples, intersect(rownames(bin), names(group_vec)))
    g <- factor(group_vec[s], levels = c(g1, g2))
    present <- as.numeric(bin[s, tx]) > 0
    prev <- as.numeric(tapply(present, g, function(z) mean(z, na.rm = TRUE))) * 100
    data.frame(
      taxa = tx,
      group = factor(c(g1, g2), levels = c(g1, g2)),
      prevalence = prev,
      signed = ifelse(c(g1, g2) == g1, prev, -prev),
      stringsAsFactors = FALSE
    )
  })
  bind_rows(out)
}

# Mean non-zero CLR difference (g1 minus g2).
# Stars are attached later from the genome-wide LM (non-zero CLR ~ EC_yesno / INR_YESNO).
clr_mean_diff <- function(taxa, samples, group_vec, g1, g2) {
  plot_samples <- intersect(samples, names(group_vec)[group_vec %in% c(g1, g2)])
  out <- lapply(taxa, function(tx) {
    empty <- data.frame(
      taxa = tx, diff_mean = NA_real_, higher_in = NA_character_,
      stringsAsFactors = FALSE
    )
    clr_mat <- clr_mat_for(tx)
    if (!tx %in% colnames(clr_mat)) return(empty)
    keep_s <- intersect(plot_samples, rownames(clr_mat))
    x <- as.numeric(clr_mat[keep_s, tx])
    g <- group_vec[keep_s]
    d1 <- x[g == g1]; d0 <- x[g == g2]
    d1 <- d1[is.finite(d1)]; d0 <- d0[is.finite(d0)]
    if (length(d1) == 0 || length(d0) == 0) return(empty)
    dm <- mean(d1) - mean(d0)
    data.frame(
      taxa = tx, diff_mean = dm,
      higher_in = ifelse(dm >= 0, g1, g2),
      stringsAsFactors = FALSE
    )
  })
  bind_rows(out)
}

ec_group_vec <- setNames(
  ifelse(as.character(pheno_EC$EC_yesno) == "EC", "EC", "NP"),
  rownames(pheno_EC)
)
ir_group_vec <- setNames(
  ifelse(as.character(pheno_IR$INR_YESNO) == "Non-Responder", "INR", "IR"),
  rownames(pheno_IR)
)

# Group LM and LRM for EC vs NP and INR vs IR.
# The grouping column is a factor. Latest CD4 count is not a covariate.
# FDR is computed across the full genus or SGB table.
pheno_EC$EC_yesno <- factor(pheno_EC$EC_yesno)
pheno_IR$INR_YESNO <- factor(pheno_IR$INR_YESNO)

run_lm_full <- function(y_mat, pheno_df, cov_df, samples, group_col) {
    lm_phenotype(
    y_mat[samples, , drop = FALSE],
    pheno_df[samples, group_col, drop = FALSE],
    cov_df[samples, , drop = FALSE],
    colnames(cov_df)
  )
}
run_lrm_full <- function(y_mat, pheno_df, cov_df, samples, group_col) {
    lrm_phenotypes(
    y_mat[samples, , drop = FALSE],
    pheno_df[samples, group_col, drop = FALSE],
    cov_df[samples, , drop = FALSE],
    colnames(cov_df)
  )
}

extract_da <- function(combined, taxa) {
  tax <- as.character(
    if ("taxa" %in% names(combined)) combined$taxa
    else if ("Taxa" %in% names(combined)) combined$Taxa
    else combined$Y
  )
  p <- if ("p" %in% names(combined)) combined$p else combined$P_value
  fdr <- if ("fdr.p" %in% names(combined)) combined$fdr.p else combined$P_adjusted_FDR
  sgb_id <- function(x) sub(".*?(SGB[0-9]+).*", "t__\\1", as.character(x), ignore.case = TRUE)
  i <- match(as.character(taxa), tax)
  miss <- is.na(i)
  if (any(miss)) i[miss] <- match(sgb_id(taxa[miss]), sgb_id(tax))
  data.frame(
    taxa = as.character(taxa),
    p = as.numeric(p[i]),
    fdr.p = as.numeric(fdr[i]),
    stringsAsFactors = FALSE
  )
}

sample_ec_lrm <- Reduce(intersect, list(
  rownames(pheno_EC), rownames(pheno_complete_ec),
  rownames(sp_re_metaphlan_genus_binary)
))
sample_ir_lrm_genus <- Reduce(intersect, list(
  rownames(pheno_IR), rownames(pheno_complete_ir),
  rownames(sp_re_metaphlan_genus_binary)
))

check_da <- function(da, must_fdr, label) {
  bad <- !is.finite(da$p) | !is.finite(da$fdr.p)
  if (any(bad)) {
    stop(label, " missing p/FDR for: ", paste(da$taxa[bad], collapse = ", "))
  }
  miss <- da$taxa %in% must_fdr & da$fdr.p >= 0.1
  if (any(miss)) {
    warning(
      label, " FDR is not < 0.1 for expected hits: ",
      paste(da$taxa[miss], collapse = ", "),
      " (p = ", paste(signif(da$p[miss], 3), collapse = ", "),
      "; FDR = ", paste(signif(da$fdr.p[miss], 3), collapse = ", "), ")."
    )
  }
}

lm_ec_gen_full <- run_lm_full(
  sp_re_metaphlan_clr_genus_non_zero, pheno_EC, pheno_complete_ec,
  sample_ec_genus, "EC_yesno"
)
da_lm_ec <- extract_da(lm_ec_gen_full, taxa_ec_lm_genus)

lm_ir_gen_full <- run_lm_full(
  sp_re_metaphlan_clr_genus_non_zero, pheno_IR, pheno_complete_ir,
  sample_ir_genus, "INR_YESNO"
)
da_lm_inr_gen <- extract_da(lm_ir_gen_full, taxa_inr_lm_genus)

lm_ir_sgb_full <- run_lm_full(
  sp_re_metaphlan_clr_non_zero, pheno_IR, pheno_complete_ir,
  sample_ir_sgb, "INR_YESNO"
)
da_lm_inr_sgb <- extract_da(lm_ir_sgb_full, c(taxa_inr_lm_sgb, taxa_inr_lrm_sgb))

lrm_ec_gen_full <- run_lrm_full(
  sp_re_metaphlan_genus_binary, pheno_EC, pheno_complete_ec,
  sample_ec_lrm, "EC_yesno"
)
lrm_ir_gen_full <- run_lrm_full(
  sp_re_metaphlan_genus_binary, pheno_IR, pheno_complete_ir,
  sample_ir_lrm_genus, "INR_YESNO"
)
lrm_ir_sgb_full <- run_lrm_full(
  sp_re_metaphlan_binary, pheno_IR, pheno_complete_ir,
  sample_ir_lrm, "INR_YESNO"
)

check_da(da_lm_ec, taxa_ec_lm_genus, "EC genus LM")
check_da(da_lm_inr_gen, taxa_inr_lm_genus, "INR genus LM")
check_da(da_lm_inr_sgb, taxa_inr_lm_sgb, "INR SGB LM")

da_lrm_ec <- extract_da(lrm_ec_gen_full, taxa_ec_lm_genus)
da_lm_inr <- rbind(da_lm_inr_gen, da_lm_inr_sgb)
da_lrm_inr <- rbind(
  extract_da(lrm_ir_gen_full, taxa_inr_lm_genus),
  extract_da(lrm_ir_sgb_full, c(taxa_inr_lm_sgb, taxa_inr_lrm_sgb))
)
check_da(da_lrm_inr, taxa_inr_lrm_sgb, "INR SGB LRM")
star_preview <- function(da, panel) {
  ann <- annotate_sig(da$p, da$fdr.p)
  cbind(panel = panel, da, glyph = ann$signif, from = ann$sig_level)
}
print(rbind(
  star_preview(da_lm_ec,  "EC"),
  star_preview(da_lm_inr, "INR")
), row.names = FALSE)
print(rbind(
  star_preview(da_lrm_ec,  "EC"),
  star_preview(da_lrm_inr, "INR")
), row.names = FALSE)

diff_ec <- clr_mean_diff(
  taxa_ec_lm_genus,
  sample_ec_genus,
  ec_group_vec, "EC", "NP"
)
diff_inr <- clr_mean_diff(
  taxa_inr_all,
  unique(c(sample_ir_genus, sample_ir_sgb, sample_ir_lrm)),
  ir_group_vec, "INR", "IR"
)

prep_diff <- function(df, taxa_keep, g_levels) {
  df$taxa_lab <- factor(unname(taxon_labs[as.character(df$taxa)]),
                        levels = unname(taxon_labs[taxa_keep]))
  df$higher_in <- factor(df$higher_in, levels = g_levels)
  df
}
diff_ec  <- prep_diff(diff_ec,  taxa_ec_lm_genus, names(cols_ec))
diff_inr <- prep_diff(diff_inr, taxa_inr_all,     names(cols_ir))

prev_ec <- taxon_prev(
  taxa_ec_lm_genus,
  rownames(pheno_EC),
  ec_group_vec, "EC", "NP"
)
prev_inr <- taxon_prev(
  taxa_inr_all,
  rownames(pheno_IR),
  ir_group_vec, "INR", "IR"
)
prev_ec$taxa_lab  <- factor(unname(taxon_labs[as.character(prev_ec$taxa)]),
                            levels = unname(taxon_labs[taxa_ec_lm_genus]))
prev_inr$taxa_lab <- factor(unname(taxon_labs[as.character(prev_inr$taxa)]),
                            levels = unname(taxon_labs[taxa_inr_all]))

attach_da_sig <- function(df, da) {
  m <- match(as.character(df$taxa), as.character(da$taxa))
  df$p <- da$p[m]
  df$fdr.p <- da$fdr.p[m]
  ann <- annotate_sig(df$p, df$fdr.p)
  df$signif <- ann$signif
  df$sig_level <- ann$sig_level
  df
}

star_layers <- function(star_df) {
  list(
    geom_text(
      data = star_df[star_df$sig_level == "FDR", ],
      aes(x = x_pos, y = star_y, label = signif),
      inherit.aes = FALSE, colour = "black", size = 4
    )
  )
}

add_x <- function(df, taxa_keep, tax_col = "taxa") {
  df$x_pos <- match(as.character(df[[tax_col]]), taxa_keep)
  df
}
prev_ec  <- add_x(attach_da_sig(prev_ec,  da_lrm_ec),  taxa_ec_lm_genus)
prev_inr <- add_x(attach_da_sig(prev_inr, da_lrm_inr), taxa_inr_all)
diff_ec  <- add_x(attach_da_sig(diff_ec,  da_lm_ec),   taxa_ec_lm_genus)
diff_inr <- add_x(attach_da_sig(diff_inr, da_lm_inr),  taxa_inr_all)

dmax <- max(abs(c(diff_ec$diff_mean, diff_inr$diff_mean)), na.rm = TRUE)
if (!is.finite(dmax) || dmax == 0) dmax <- 1
star_pad <- 0.12 * dmax
diff_lim <- range(c(diff_ec$diff_mean, diff_inr$diff_mean), na.rm = TRUE)
diff_lim <- c(min(diff_lim[1], 0) - 1.35 * star_pad,
              max(diff_lim[2], 0) + 1.35 * star_pad)

prev_star_df <- function(df) {
  d <- df %>%
    dplyr::group_by(.data$taxa, .data$x_pos, .data$signif, .data$sig_level) %>%
    dplyr::summarise(
      star_y = pmin(108, suppressWarnings(max(.data$signed, na.rm = TRUE)) + 8),
      .groups = "drop"
    )
  d[d$sig_level != "None" & is.finite(d$star_y), ]
}
diff_star_df <- function(df) {
  d <- df
  d$star_y <- ifelse(
    is.na(d$diff_mean), NA_real_,
    d$diff_mean + ifelse(d$diff_mean >= 0, 1, -1) * star_pad
  )
  d[d$sig_level != "None" & is.finite(d$star_y),
    c("x_pos", "star_y", "signif", "sig_level")]
}

stars_prev_ec  <- prev_star_df(prev_ec)
stars_prev_inr <- prev_star_df(prev_inr)
stars_diff_ec  <- diff_star_df(diff_ec)
stars_diff_inr <- diff_star_df(diff_inr)

x_scale <- function(taxa_keep, show_labels = FALSE) {
  labs <- unname(taxon_labs[taxa_keep])
  scale_x_continuous(
    limits = c(0.5, length(taxa_keep) + 0.5),
    breaks = seq_along(taxa_keep),
    labels = if (show_labels) labs else NULL,
    expand = c(0, 0)
  )
}

make_prev <- function(df, taxa_keep, cols, title, subtitle, show_y = TRUE,
                      lrm_line = FALSE, star_df = NULL) {
  p <- ggplot(df, aes(x = x_pos, y = signed, fill = group)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.35)
  if (!is.null(star_df) && nrow(star_df) > 0) p <- p + star_layers(star_df)
  p <- p +
    scale_fill_manual(values = cols, name = NULL, drop = FALSE) +
    x_scale(taxa_keep, show_labels = FALSE) +
    scale_y_continuous(
      breaks = seq(-100, 100, 50),
      labels = function(v) abs(v)
    ) +
    coord_cartesian(ylim = c(-100, 115), expand = FALSE) +
    labs(x = NULL, y = "Prevalence (%)",
         title = title, subtitle = subtitle) +
    theme_minimal() +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey30"),
      axis.text.x   = element_blank(),
      axis.ticks.x  = element_blank(),
      axis.text.y   = element_text(size = 9),
      panel.grid.major.x = element_blank(),
      plot.margin   = margin(4, 4, 2, 4)
    )
  if (lrm_line && any(taxa_keep == "t__SGB14754")) {
    p <- p + geom_vline(
      xintercept = sum(taxa_keep != "t__SGB14754") + 0.5,
      colour = "grey55", linewidth = 0.4
    )
  }
  if (!show_y) p <- p + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  p
}

make_diff <- function(df, taxa_keep, cols, show_y = TRUE, lrm_line = FALSE,
                      star_df = NULL) {
  p <- ggplot(df, aes(x = x_pos, y = diff_mean, colour = higher_in)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.35) +
    geom_segment(aes(xend = x_pos, y = 0, yend = diff_mean),
                 linewidth = 0.55, na.rm = TRUE) +
    geom_point(size = 2.6, na.rm = TRUE)
  if (!is.null(star_df) && nrow(star_df) > 0) p <- p + star_layers(star_df)
  p <- p +
    scale_colour_manual(values = cols, name = "Higher CLR in", drop = FALSE) +
    x_scale(taxa_keep, show_labels = FALSE) +
    coord_cartesian(ylim = diff_lim) +
    labs(x = NULL, y = "Δ CLR") +
    theme_minimal() +
    theme(
      axis.text.x   = element_blank(),
      axis.ticks.x  = element_blank(),
      axis.text.y   = element_text(size = 9),
      panel.grid.major.x = element_blank(),
      plot.margin   = margin(2, 4, 2, 4)
    )
  if (lrm_line && any(taxa_keep == "t__SGB14754")) {
    p <- p + geom_vline(
      xintercept = sum(taxa_keep != "t__SGB14754") + 0.5,
      colour = "grey55", linewidth = 0.4
    )
  }
  if (!show_y) p <- p + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  p
}

make_heat <- function(df, taxa_keep, lrm_line = FALSE, show_y = TRUE) {
  d <- df[as.character(df$Taxa) %in% taxa_keep, ]
  d <- add_x(d, taxa_keep, tax_col = "Taxa")
  x_cols <- unname(axis_cols[taxa_keep])
  p <- ggplot(d, aes(x = x_pos, y = Phenotype, fill = Beta)) +
    geom_tile(width = 1, height = 1, color = "white") +
    geom_text(aes(label = signif, color = sig_level), size = 5) +
    fill_scale + sig_scale +
    x_scale(taxa_keep, show_labels = TRUE) +
    scale_y_discrete(limits = rev(drug_vars),
                     labels = pheno_labels[rev(drug_vars)],
                     expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x   = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   colour = x_cols, size = 9, lineheight = 0.9),
      axis.text.y   = element_text(size = 9),
      panel.grid    = element_blank(),
      plot.margin   = margin(2, 4, 4, 4)
    )
  if (lrm_line && any(taxa_keep == "t__SGB14754")) {
    p <- p + geom_vline(
      xintercept = sum(taxa_keep != "t__SGB14754") + 0.5,
      colour = "grey55", linewidth = 0.4
    )
  }
  if (!show_y) p <- p + theme(axis.text.y = element_blank())
  p
}

p_prev_ec <- make_prev(
  prev_ec, taxa_ec_lm_genus, cols_ec,
  title = "EC vs NP", subtitle = n_lab(length(sample_ec_genus)),
  star_df = stars_prev_ec
)
p_prev_inr <- make_prev(
  prev_inr, taxa_inr_all, cols_ir,
  title = "INR vs IR", subtitle = n_lab(length(sample_ir_genus)),
  show_y = FALSE, lrm_line = TRUE, star_df = stars_prev_inr
)
p_diff_ec <- make_diff(diff_ec, taxa_ec_lm_genus, cols_ec, star_df = stars_diff_ec)
p_diff_inr <- make_diff(diff_inr, taxa_inr_all, cols_ir, show_y = FALSE,
                       lrm_line = TRUE, star_df = stars_diff_inr)
p_a_ec <- make_heat(heat_within, taxa_ec_lm_genus)
p_a_inr <- make_heat(heat_within, taxa_inr_all, lrm_line = TRUE, show_y = FALSE)
p_b <- make_heat(heat_ir, taxa_order, lrm_line = TRUE) +
  labs(title = "IR", subtitle = n_lab(length(sample_genus))) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey30")
  )

# A: prevalence + Δ CLR + within-cohort heatmaps; B: Responders
# Same numeric x in each column (no axes="collect") so dots sit on heatmap tiles.
p_a <- p_prev_ec + p_prev_inr + p_diff_ec + p_diff_inr + p_a_ec + p_a_inr +
  plot_layout(
    design = "AB\nCD\nEF",
    widths = c(length(taxa_ec_lm_genus), length(taxa_inr_all)),
    heights = c(0.65, 0.65, 1.5),
    guides = "collect"
  )
fig_linked_pheno <- p_a / p_b +
  plot_layout(heights = c(2.8, 1.5), guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    caption = paste(
      "Prevalence: % of samples with the taxon present (upper bar = EC or INR, lower = NP or IR).",
      "Δ CLR: mean non-zero CLR in EC minus NP, or INR minus IR.",
      "Prevalence stars: LRM of presence vs EC vs NP or INR vs IR (genome-wide FDR).",
      "Δ CLR stars: LM of non-zero CLR vs the same grouping (genome-wide FDR).",
      "One black star if FDR < 0.1 (heatmap, prevalence and delta CLR).",
      "LM = non-zero CLR; LRM = presence/absence."
    )
  ) &
  theme(legend.position = "right")
fig_linked_pheno


## Project: Gut Microbiome in Spontaneous HIV-1 Controllers and Immunological Non Responders: a metagenomics study in 1,559 people with HIV
## Author: Ángela del Castillo Izquierdo
## Last time edited: 2026-10-07
##
## 0. Prepare working files ----------------------------------------------------
## Builds the phenotype subsets and the abundance tables used by scripts 1-4.
##
## Input (1,559 participants with complete phenotypic data):
##   sgb_relative_abundance.tsv
##   genus_relative_abundance.tsv
##   pathway_relative_abundance.tsv
##   taxonomy_sgb.tsv
##   samples_metadata.tsv
##
## The prevalence filter (>= 15 samples) was applied on all sequenced MetaPhlAn
## profiles (1,691). Analysis is then limited to these 1,559 participants, 
## the ones with complete phenotypic data.
## The full-cohort detection counts are n_detected_sequenced and
## genus/order/phylum_n_detected_sequenced in taxonomy_sgb.tsv.
## A genus is uncharacterised when its id starts with g__GGB; an SGB is
## uncharacterised when the species name starts with GGB.
## 
##
## CLR:
##   SGB and pathways: log(x + min_positive/2), centered within each sample.
##   Genus, order and phylum: microbiome::transform(..., "clr") on the
##   unfiltered aggregated table, then the prevalence filter.
##
## Packages:
##   dplyr       1.2.0
##   microbiome  1.32.0

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(microbiome))

root <- path.expand("~/Documents/Projects/HIV Elite Controllers/Files_Yue_2025_10_06/code_public")
data_dir <- file.path(root, "data")
PREV_MIN_N <- 15

read_abundance <- function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(x) <- x$SampleID
  x$SampleID <- NULL
  x[] <- lapply(x, as.numeric)
  x
}

sgb_rel <- read_abundance(file.path(data_dir, "sgb_relative_abundance.tsv"))
genus_rel <- read_abundance(file.path(data_dir, "genus_relative_abundance.tsv"))
path_rel <- read_abundance(file.path(data_dir, "pathway_relative_abundance.tsv"))
taxonomy_table <- read.delim(file.path(data_dir, "taxonomy_sgb.tsv"),
                             check.names = FALSE, stringsAsFactors = FALSE)
pheno <- read.delim(file.path(data_dir, "samples_metadata.tsv"),
                    check.names = FALSE, stringsAsFactors = FALSE)
rownames(pheno) <- pheno$SampleID

# --- Phenotype subsets --------------------------------------------------------
# pheno            all analysed participants                         N = 1559
# pheno_EC         normal progressors + persistent non-viremic EC    N = 1483
# pheno_IR         responders and non-responders, progressors only   N = 1238
# pheno_persistent HIC-persistent vs normal progressors              N = 1513
# pheno_transient  HIC-transient vs normal progressors               N = 1510
pheno_EC <- pheno[pheno$EC_combined %in% c("non_EC", "EC_persistent_nonviremic_noART"), ]
pheno_IR <- pheno[pheno$INR_YESNO %in% c("Responder", "Non-Responder") &
                    pheno$EC_yesno == "non_EC", ]
pheno_persistent <- pheno[pheno$EC_lost %in% c("non_EC", "EC_persistent"), ]
pheno_transient <- pheno[pheno$EC_lost %in% c("non_EC", "EC_transient"), ]
rownames(pheno_EC) <- pheno_EC$SampleID
rownames(pheno_IR) <- pheno_IR$SampleID
rownames(pheno_persistent) <- pheno_persistent$SampleID
rownames(pheno_transient) <- pheno_transient$SampleID

stopifnot(nrow(pheno) == 1559)
stopifnot(nrow(pheno_EC) == 1483)
stopifnot(nrow(pheno_IR) == 1238)
stopifnot(nrow(pheno_persistent) == 1513)
stopifnot(nrow(pheno_transient) == 1510)

# Centered log-ratio with pseudocount = half the smallest positive value.
# Computed on the full feature table so the geometric mean includes rare taxa,
# then columns are subset to the prevalence filter.
clr_half_minimum <- function(abundance) {
  abundance <- as.matrix(abundance)
  positive <- abundance[is.finite(abundance) & abundance > 0]
  pseudo <- min(positive) / 2
  transformed <- log(abundance + pseudo)
  transformed <- transformed - rowMeans(transformed)
  list(clr = as.data.frame(transformed, check.names = FALSE), pseudo = pseudo)
}

finish_tables <- function(rel_all, clr_all, keep) {
  keep <- intersect(keep, colnames(rel_all))
  rel <- rel_all[, keep, drop = FALSE]
  clr <- clr_all[, keep, drop = FALSE]
  clr_nz <- clr
  clr_nz[rel == 0] <- NA
  binary <- rel
  binary[binary > 0] <- 1
  list(rel = rel, binary = binary, clr = clr, clr_nz = clr_nz)
}

# --- SGB ----------------------------------------------------------------------
sgb_n <- taxonomy_table$n_detected_sequenced[match(colnames(sgb_rel), taxonomy_table$ID)]
keep_sgb <- colnames(sgb_rel)[!is.na(sgb_n) & sgb_n >= PREV_MIN_N]
sgb_clr <- clr_half_minimum(sgb_rel)
sgb_tables <- finish_tables(sgb_rel, sgb_clr$clr, keep_sgb)
sp_re_metaphlan <- sgb_tables$rel
sp_re_metaphlan_binary <- sgb_tables$binary
sp_re_metaphlan_clr <- sgb_tables$clr
sp_re_metaphlan_clr_non_zero <- sgb_tables$clr_nz

# --- Genus, order and phylum --------------------------------------------------
# Aggregated from every SGB (rare SGBs included). CLR is computed on that full
# table, then taxa detected in >= PREV_MIN_N sequenced samples are kept.
# Genus characterisation elsewhere is startsWith(id, "g__GGB").
genus_ids_unfiltered <- colnames(genus_rel)
genus_n <- taxonomy_table$genus_n_detected_sequenced[
  match(colnames(genus_rel), taxonomy_table$genus_id)
]
keep_genus <- colnames(genus_rel)[!is.na(genus_n) & genus_n >= PREV_MIN_N]
genus_clr_all <- t(microbiome::transform(t(as.matrix(genus_rel)), "clr"))
genus_clr_all <- as.data.frame(genus_clr_all, check.names = FALSE)
genus_tables <- finish_tables(genus_rel, genus_clr_all, keep_genus)
sp_re_metaphlan_genus <- genus_tables$rel
sp_re_metaphlan_genus_binary <- genus_tables$binary
sp_re_metaphlan_clr_genus <- genus_tables$clr
sp_re_metaphlan_clr_genus_non_zero <- genus_tables$clr_nz

aggregate_rank <- function(rank_name) {
  vec <- taxonomy_table[[rank_name]][match(colnames(sgb_rel), taxonomy_table$ID)]
  vec[is.na(vec) | vec == ""] <- "Unknown"
  agg <- t(rowsum(t(as.matrix(sgb_rel)), group = vec, reorder = TRUE))
  storage.mode(agg) <- "double"
  as.data.frame(agg, check.names = FALSE)
}
filter_rank <- function(rel, n_by_taxon) {
  clr_all <- t(microbiome::transform(t(as.matrix(rel)), "clr"))
  clr_all <- as.data.frame(clr_all, check.names = FALSE)
  n <- n_by_taxon[colnames(rel)]
  keep <- colnames(rel)[!is.na(n) & n >= PREV_MIN_N]
  finish_tables(rel, clr_all, keep)
}
order_tables <- filter_rank(
  aggregate_rank("Order"),
  setNames(taxonomy_table$order_n_detected_sequenced, taxonomy_table$Order)
)
phylum_tables <- filter_rank(
  aggregate_rank("Phylum"),
  setNames(taxonomy_table$phylum_n_detected_sequenced, taxonomy_table$Phylum)
)
sp_re_metaphlan_order <- order_tables$rel
sp_re_metaphlan_order_binary <- order_tables$binary
sp_re_metaphlan_clr_order <- order_tables$clr
sp_re_metaphlan_clr_order_non_zero <- order_tables$clr_nz
sp_re_metaphlan_phylum <- phylum_tables$rel
sp_re_metaphlan_phylum_binary <- phylum_tables$binary
sp_re_metaphlan_clr_phylum <- phylum_tables$clr
sp_re_metaphlan_clr_phylum_non_zero <- phylum_tables$clr_nz

# --- Pathways (no prevalence filter) ------------------------------------------
path_clr <- clr_half_minimum(path_rel)
path_tables <- finish_tables(path_rel, path_clr$clr, colnames(path_rel))
re_path <- path_tables$rel
re_path_binary <- path_tables$binary
re_path_clr <- path_tables$clr
re_path_clr_non_zero <- path_tables$clr_nz

cat("SGB pseudocount:", sgb_clr$pseudo,
    "| kept", ncol(sp_re_metaphlan), "of", ncol(sgb_rel),
    "(filter counted on all sequenced samples; analysis N =", nrow(pheno), ")\n")
cat("Genus kept", ncol(sp_re_metaphlan_genus), "of", ncol(genus_rel), "\n")
cat("Order kept", ncol(sp_re_metaphlan_order),
    "| Phylum kept", ncol(sp_re_metaphlan_phylum), "\n")
cat("Pathway pseudocount:", path_clr$pseudo,
    "| pathways", ncol(re_path), "(no prevalence filter)",
    "| samples", nrow(re_path), "\n")
stopifnot(ncol(sp_re_metaphlan) == 1599)
stopifnot(ncol(sp_re_metaphlan_genus) == 986)
stopifnot(ncol(re_path) == ncol(path_rel))
stopifnot(abs(sgb_clr$pseudo - 5e-8) < 1e-12)
cat("EC vs NP:\n"); print(table(pheno_EC$EC_combined))
cat("INR vs IR (progressors only):\n"); print(table(pheno_IR$INR_YESNO))
cat("HIC-persistent vs NP:\n"); print(table(pheno_persistent$EC_lost))
cat("HIC-transient vs NP:\n"); print(table(pheno_transient$EC_lost))

save(pheno, pheno_EC, pheno_IR, pheno_persistent, pheno_transient,
     sp_re_metaphlan, sp_re_metaphlan_binary,
     sp_re_metaphlan_clr, sp_re_metaphlan_clr_non_zero,
     sp_re_metaphlan_genus, sp_re_metaphlan_genus_binary,
     sp_re_metaphlan_clr_genus, sp_re_metaphlan_clr_genus_non_zero,
     sp_re_metaphlan_order, sp_re_metaphlan_order_binary,
     sp_re_metaphlan_clr_order, sp_re_metaphlan_clr_order_non_zero,
     sp_re_metaphlan_phylum, sp_re_metaphlan_phylum_binary,
     sp_re_metaphlan_clr_phylum, sp_re_metaphlan_clr_phylum_non_zero,
     re_path, re_path_binary, re_path_clr, re_path_clr_non_zero,
     taxonomy_table, genus_ids_unfiltered, PREV_MIN_N,
     file = file.path(data_dir, "working_objects.RData"))
cat("Saved", file.path(data_dir, "working_objects.RData"), "\n")

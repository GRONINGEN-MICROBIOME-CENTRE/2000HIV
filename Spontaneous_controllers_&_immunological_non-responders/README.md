# Gut microbiome in spontaneous HIV-1 controllers and immunological non-responders: a metagenomics study in 1,559 people with HIV
Code accompanying the publication:

*Del Castillo-Izquierdo, Á., Zhang, Y., et al. (2026)* <br>
*Gut microbiome in spontaneous HIV-1 controllers and immunological non-responders: a metagenomics study in 1,559 people with HIV.*

## Contents
Run the scripts in the order below. Script 0 prepares the tables used by scripts 1–4.  <br>
Shotgun metagenomic profiles of faecal samples from the 2000HIV cohort, and the phenotype data for the same participants, are deposited separately.

0. [Prepare working files](0_Prepare_working_files.R) <br>
   Builds the phenotype subsets and the abundance tables at SGB, genus, order, phylum and pathway level: relative abundance, presence/absence, and centered log-ratio.

   [Association functions](R_functions.R)  <br>
   Inverse-normal transformation, linear models and logistic models used by scripts 1–4.

2. [Microbial profiling and diversity](1_Microbial_profiling_diversity.R) <br>
   Characterised and uncharacterised SGBs and genera, before and after the prevalence filter. PERMANOVA and alpha diversity for elite controllers versus non-progressors, immunological non-responders versus responders, and persistent or transient controllers versus non-progressors. Genus-level principal component analysis for the first two contrasts.

3. [Differential abundance](2_Differencial_abundance_analysis.R) <br>
   Linear models of non-zero centered log-ratio abundance, and logistic models of presence/absence, for the same four contrasts, at SGB, genus and pathway level. Elite controllers versus non-progressors are also fitted with latest CD4 count added to the covariates, at SGB and genus level only.
   
4. [Link taxa to phenotypes](3_Link_taxa_phenotypes.R) <br>
   Associates the taxa that were significant in the elite-controller and non-responder comparisons with clinical phenotypes.

5. [Extra analysis](4_Extra_analysis.R) <br>
   Prevotella/Bacteroides and Bacteroidales/Clostridiales abundance ratios, and a phylum-level comparison of elite controllers versus non-progressors.


# Load library
library(GCIM.GWEIS.Z)
# Setup paths to external software and LDSC reference files
plink <- "<path>/plink2"
hm3_snps <- "<path>/w_hm3.snplist"
ld_scores <- "<path>/eur_w_ld_chr/"

# Step 1: Run GWAS for exposure in the PRS training sample
# This estimates SNP main effects for the exposure trait.
# The resulting GWAS summary statistics will be used as weights to construct the PRS.
cat("Step 1: Running discovery GWAS...\n")
gwas_res <- b_gwas(
  plink_path = plink,
  dis_mydata = "mydata.PRStrained",
  dis_cov_file = "trait1_PRStrained_cov_b.txt",
  threads = 40
)

# Step 2: Compute PRS in the analysis sample
# The GWAS results from Step 1 are used to calculate the exposure PRS in the independent analysis sample.
cat("Step 2: Computing PRS in target...\n")
prs_res <- prs_scores(
  plink_path = plink,
  tar_mydata = "mydata.analysis",
  score_file = gwas_res,
  threads = 40
)

# Step 3: Replace the observed exposure with the PRS
# In GCIM-GWEIS, the observed exposure is replaced by its genetically predicted value (PRS).
# This helps reduce bias from reverse causal direction.
cat("Step 3: Replacing exposure with PRS...\n")
replaced_res <- replace_covariate_with_prs(
  dis_cov_file = "trait1_analysis_cov_b.txt",
  prs_file = prs_res,
  on_missing = "stop"
)

# Step 4: Run GCIM-GWEIS
# This tests SNP-by-PRS interaction effects on the outcome trait.
# The PRS is used as the interaction covariate instead of the observed exposure.
cat("Step 4: Running GWEIS...\n")
gweis_res <- q_gweis(
  plink_path = plink,
  tar_mydata = "mydata.analysis",
  tar_pheno_file = "trait2_analysis_out.txt",
  tar_covar_file = replaced_res,
  int_covar_index = 1,
  threads = 40
)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# Save the GCIM-GWEIS output file(optional)
glm_df <- read.table(gweis_res$glm_file, header = TRUE, comment.char = "", stringsAsFactors = FALSE)
write.table(glm_df, file = "trait2_out_q_gcim_gweis.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

# Step 5: Munge GCIM-GWEIS summary statistics for LDSC
# This prepares the GCIM-GWEIS summary statistics for LDSC analysis.
# The HapMap3 SNP list is used to restrict the analysis to high-quality SNPs.
cat("Step 5: Munging for LDSC...\n")
munge_res <- munge_ldsc_gcim(
  munge_path  = "<path>/ldsc/munge_sumstats.py",
  glm_file    = gweis_res,
  hm3_snplist = hm3_snps,
  python      = "<path>/anaconda3/envs/ldsc/bin/python"
)

# Step 6: Estimate the LDSC intercept
# LDSC is used to estimate the intercept of the GCIM-GWEIS test statistics.
# This intercept captures inflation in the genome-wide test statistics.
cat("Step 6: Computing LDSC intercept...\n")
ldsc_res <- ldsc_h2_gcim(
  ldsc_path = "<path>/ldsc/ldsc.py",
  python = "<path>/anaconda3/envs/ldsc/bin/python",
  munged_sumstats = munge_res,
  ref_ld_chr = ld_scores
)
# Step 7: Adjust GCIM-GWEIS Z-scores
# The GCIM-GWEIS Z-statistics are divided by the square root of the LDSC intercept.
# This produces the final GCIM-GWEIS-Z results with corrected test statistics.
# Step 7: Adjust Z-scores
cat("Step 7: Adjusting Z-scores...\n")
final_res <- gcim_z_adjust(
  glm_file = gweis_res,
  intercept_file = ldsc_res
)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# fully adjusted outputs 
file.copy(final_res$output_file, "trait2_out_q_gcim_gweis-z.txt", overwrite = TRUE)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
# Report final results
cat("\n===== PIPELINE COMPLETE =====\n")
cat("Total SNPs analyzed:", final_res$n_snps, "\n")
cat("LDSC intercept:", sprintf("%.4f", final_res$intercept_used), "\n")
cat("Correction factor:", sprintf("%.4f", final_res$correction_factor), "\n")
cat("Significant SNPs (p<0.05):\n")
cat("  Before adjustment:", final_res$n_sig_original, "\n")
cat("  After adjustment:", final_res$n_sig_adjusted, "\n")
cat("\nFinal results saved to:", final_res$output_file, "\n")

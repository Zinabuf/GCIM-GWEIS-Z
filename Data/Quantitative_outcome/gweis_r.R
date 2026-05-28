# Load library
library(GCIM.GWEIS.Z)
# Setup paths to external software and LDSC reference files
plink <- "<plink_path>/plink2"
hm3_snps <- "<path>/w_hm3.snplist"
ld_scores <- "<path>/eur_w_ld_chr/"

# Step 4: Run conventioan GWEIS for the reverse direction
# This tests SNP-by-environmental exposure interaction effects on the outcome trait.

# Step 4: GWEIS
cat("Step 4: Running GWEIS...\n")
gweis_res <- q_gweis(
  plink_path = plink,
  tar_mydata = "mydata.analysis",
  tar_pheno_file = "trait2_analysis_out.txt",
  tar_covar_file = "trait1_analysis_cov.txt",
  int_covar_index = 1,
  threads = 40
)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# Save the GCIM-GWEIS output file(optional)
glm_df <- read.table(gweis_res$glm_file, header = TRUE, comment.char = "", stringsAsFactors = FALSE)
write.table(glm_df, file = "trait2_out_q_gweis.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

# Step 5: Munge convetional GWEIS summary statistics for LDSC
# This prepares the conventional GWEIS summary statistics for LDSC analysis.
# The HapMap3 SNP list is used to restrict the analysis to high-quality SNPs.
cat("Step 5: Munging for LDSC...\n")
munge_res <- munge_ldsc_gcim(
  munge_path  = "/data/alh-admzw/ldsc/munge_sumstats.py",
  glm_file    = gweis_res,
  hm3_snplist = hm3_snps,
  python      = "/data/alh-admzw/anaconda3/envs/ldsc/bin/python"
)
# Step 6: Estimate the LDSC intercept
# LDSC is used to estimate the intercept of the conventional GWEIS test statistics.
# This intercept captures inflation in the genome-wide test statistics.
cat("Step 6: Computing LDSC intercept...\n")
ldsc_res <- ldsc_h2_gcim(
  ldsc_path = "/data/alh-admzw/ldsc/ldsc.py",
  python = "/data/alh-admzw/anaconda3/envs/ldsc/bin/python",
  munged_sumstats = munge_res,
  ref_ld_chr = ld_scores
)

# Step 7: Adjust conventional GWEIS Z-scores
# The conventional GWEIS Z-statistics are divided by the square root of the LDSC intercept.
# This produces the final GWEIS-Z results with corrected test statistics for the reverse causal directions.
cat("Step 7: Adjusting Z-scores...\n")
final_res <- gcim_z_adjust(
  glm_file = gweis_res,
  intercept_file = ldsc_res
)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# Save the final GWEIS-Z adjusted output 
file.copy(final_res$output_file, "trait2_out_q_gweis-z.txt", overwrite = TRUE)
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

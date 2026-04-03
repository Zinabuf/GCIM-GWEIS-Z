library(GCIM.GWEIS.Z)
# Setup paths
# Setup paths
plink <- "<path>/plink2"
munge_path <- "<path>/munge_sumstats.py"
python     <- "<path>/anaconda3/envs/ldsc/bin/python"

#munge_path <- "<path>/ldsc/munge_sumstats.py"
ldsc_path <- "<path>/ldsc/ldsc.py"
hm3_snps <- "<path>/w_hm3.snplist"
ld_scores <- "<path>/eur_w_ld_chr/"
python <- "<path>/anaconda3"

# Step 4: GWEIS
cat("Step 4: Running GWEIS...\n")
gweis_res <- q_gweis(
  plink_path = plink,
  tar_mydata = "mydata.analysis",
  tar_pheno_file = "trait2_analysis_out.txt",
  tar_covar_file = "trait1_analysis_cov_b.txt",
  int_covar_index = 1,
  threads = 40
)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
glm_df <- read.table(gweis_res$glm_file, header = TRUE, comment.char = "", stringsAsFactors = FALSE)
write.table(glm_df, file = "trait2_out_q_gweis.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

#write.table(gweis_res, "gcimt1.PHENO1.glm.linear", quote = F, row.names = F, sep = " ")
# Step 5: Munge for LDSC
cat("Step 5: Munging for LDSC...\n")
munge_res <- munge_ldsc_gcim(
  munge_path  = "<path>/ldsc/munge_sumstats.py",
  glm_file    = gweis_res,
  hm3_snplist = hm3_snps,
  python      = "<path>/anaconda3/envs/ldsc/bin/python"
)
# Step 6: LDSC heritability
cat("Step 6: Computing LDSC intercept...\n")
ldsc_res <- ldsc_h2_gcim(
  ldsc_path = "<path>/ldsc/ldsc.py",
  python = "<path>/anaconda3/envs/ldsc/bin/python",
  munged_sumstats = munge_res,
  ref_ld_chr = ld_scores
)

# Step 7: Adjust Z-scores
cat("Step 7: Adjusting Z-scores...\n")
final_res <- gcim_z_adjust(
  glm_file = gweis_res,
  intercept_file = ldsc_res
)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# fully adjusted outputs 
zdf <- read.table(final_res$output_file, header = TRUE, stringsAsFactors = FALSE)
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

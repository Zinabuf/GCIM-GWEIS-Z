library(GCIM.GWEIS.Z)
# Setup paths
plink <- "/data/alh-admzw/bscaling/basic_sim/plink2"
munge_path <- "/data/alh-admzw/ldsc/munge_sumstats.py"
python     <- "/data/alh-admzw/anaconda3/envs/ldsc/bin/python"

#munge_path <- "/home/567/zw6700/ldsc/munge_sumstats.py"
ldsc_path <- "/data/alh-admzw/ldsc/ldsc.py"
hm3_snps <- "/data/alh-admzw/ldsc_chapter3/real_data/w_hm3.snplist"
ld_scores <- "/data/alh-admzw/ldsc_chapter3/real_data/eur_w_ld_chr/"
python <- "/data/alh-admzw/anaconda3"
# Discovery data
dis_geno <- "/data/alh-admzw/ldsc_chapter3/real_data/large_target/revised/true/bil_bmi/package/toydata/mydata.dis"
#dis_pheno <- "/scratch/ny48/zw6700/Chpater_3/null.corr/scaled/null_corr_rxe05/g08/g8r8/phend"

# Target data
tar_geno <- "/data/alh-admzw/ldsc_chapter3/real_data/large_target/revised/true/bil_bmi/package/toydata/mydata.tar"
#tar_pheno <- "/scratch/ny48/zw6700/Chpater_3/null.corr/scaled/null_corr_rxe05/g08/g8r8/phent"
#tar_covar <- "/scratch/ny48/zw6700/Chpater_3/null.corr/scaled/null_corr_rxe05/g08/g8r8/covt"

# Step 1: Discovery GWAS
cat("Step 1: Running discovery GWAS...\n")
gwas_res <- q_gwas(
  plink_path = plink,
  dis_mydata = "mydata.PRStrained",
  dis_cov_file = "trait2_PRStrained_cov.txt",
  threads = 40
)

# Step 2: Compute PRS
cat("Step 2: Computing PRS in target...\n")
prs_res <- prs_scores(
  plink_path = plink,
  tar_mydata = "mydata.analysis",
  score_file = gwas_res,
  threads = 40
)

# Step 3: Replace exposure with PRS
cat("Step 3: Replacing exposure with PRS...\n")
replaced_res <- replace_covariate_with_prs(
  dis_cov_file = "trait2_analysis_cov.txt",
  prs_file = prs_res,
  on_missing = "stop"
)

# Step 4: GWEIS
cat("Step 4: Running GWEIS...\n")
gweis_res <- q_gweis(
  plink_path = plink,
  tar_mydata = "mydata.analysis",
  tar_pheno_file = "trait1_analysis_out.txt",
  tar_covar_file = replaced_res,
  int_covar_index = 1,
  threads = 40
)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
glm_df <- read.table(gweis_res$glm_file, header = TRUE, comment.char = "", stringsAsFactors = FALSE)
write.table(glm_df, file = "trait1_out_q_gcim_gweis.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

#write.table(gweis_res, "gcimt1.PHENO1.glm.linear", quote = F, row.names = F, sep = " ")
# Step 5: Munge for LDSC
cat("Step 5: Munging for LDSC...\n")
munge_res <- munge_ldsc_gcim(
  munge_path  = "/data/alh-admzw/ldsc/munge_sumstats.py",
  glm_file    = gweis_res,
  hm3_snplist = hm3_snps,
  python      = "/data/alh-admzw/anaconda3/envs/ldsc/bin/python"
)
# Step 6: LDSC heritability
cat("Step 6: Computing LDSC intercept...\n")
ldsc_res <- ldsc_h2_gcim(
  ldsc_path = "/data/alh-admzw/ldsc/ldsc.py",
  python = "/data/alh-admzw/anaconda3/envs/ldsc/bin/python",
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
file.copy(final_res$output_file, "trait1_out_q_gcim_gweis-z.txt", overwrite = TRUE)
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

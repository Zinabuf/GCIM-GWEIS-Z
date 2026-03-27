library(GCIM.GWEIS.Z)
# Setup paths
plink <- "<path>/plink2"
#munge_path <- "<path>/ldsc/munge_sumstats.py"
#python     <- "<path>/anaconda3/envs/ldsc/bin/python"

#munge_path <- "<path>/ldsc/munge_sumstats.py"
#ldsc_path <- "<path>/ldsc/ldsc.py"
hm3_snps <- "<path>/ldsc_chapter3/real_data/w_hm3.snplist"
ld_scores <- "<path>/ldsc_chapter3/real_data/eur_w_ld_chr/"
#python <- "<path>/anaconda3"

# Step 1: Discovery GWAS
cat("Step 1: Running discovery GWAS...\n")
gwas_res <- q_gwas(
  plink_path = plink,
  dis_mydata = "mydata.dis",
  dis_cov_file = "bmi_dis_cov.txt",
  threads = 40
)

# Step 2: Compute PRS
cat("Step 2: Computing PRS in target...\n")
prs_res <- prs_scores(
  plink_path = plink,
  tar_mydata = "mydata.tar",
  score_file = gwas_res,
  threads = 40
)

# Step 3: Replace exposure with PRS
cat("Step 3: Replacing exposure with PRS...\n")
replaced_res <- replace_covariate_with_prs(
  dis_cov_file = "bmi_tar_cov.txt",
  prs_file = prs_res,
  on_missing = "stop"
)

# Step 4: GWEIS
cat("Step 4: Running GWEIS...\n")
gweis_res <- q_gweis(
  plink_path = plink,
  tar_mydata = "mydata.tar",
  tar_pheno_file = "tBil_tar_out.txt",
  tar_covar_file = replaced_res,
  int_covar_index = 1,
  threads = 40
)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
glm_df <- read.table(gweis_res$glm_file, header = TRUE, comment.char = "", stringsAsFactors = FALSE)
write.table(glm_df, file = "bmi_tBil_gcim-gweis.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
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
file.copy(final_res$output_file, "bmi_tBil_gcim-gweis-z.txt", overwrite = TRUE)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
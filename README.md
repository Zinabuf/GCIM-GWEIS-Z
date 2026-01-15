---
# <h1 align="center">GCIM-GWEIS-Z</h1>
---

The Genetic Causality Inference Model (GCIM) for genome-wide-by-environment interaction studies (GWEIS), incorporating adjusted z-scores (Z), is a statistical framework designed to infer the causal direction of SNP-by-environment interactions and to correct inflation in GxE effects arising from heteroscedasticity.
#### Authors: Zinabu Fentaw and S.Hong Lee

   
## 1. Package installation 
From GitHub 

~~~
library(devtools)
install_github("Zinabuf/GCIM-GWEIS-Z")
~~~

## 2. Load the library


~~~
library(GCIM.GWEIS.Z)
~~~

## 3. Analysis pipeline
### 1. Quantitative Outcomes with Quantitative Exposures 

~~~
#load library
library(GCIM.GWEIS.Z)
# Set plink path
plink_path <- "<plink_path>/plink2"
# Set munge_sumstats.py path
menge_path <- "<menge_path>/munge_sumstats.py"
# Set ldsc.py path
ldsc_path <- "<ldsc_path>/ldsc.py"
hm3_snps <- "<snplst_path>/w_hm3.snplist"
ld_scores <- "<ld_ref_path>/eur_w_ld_chr/"

# Step 1: Discovery GWAS
cat("Step 1: Running discovery GWAS...\n")
gwas_res <- q_gwas(plink_path = plink, "dis_mydata", dis_cov_file = "dis_cov", threads = 40)
# Step 2: Compute PRS
cat("Step 2: Computing PRS in target...\n")
prs_res <- prs_scores(plink_path = plink, "tar_mydata", score_file = gwas_res, threads = 40)

# Step 3: Replace exposure with PRS
cat("Step 3: Replacing exposure with PRS...\n")
replaced_res <- replace_covariate_with_prs("tar_covar", prs_file = prs_res, on_missing = "stop")
# Step 4: GWEIS
cat("Step 4: Running GWEIS...\n")
gweis_res <- q_gweis(plink_path = plink, "tar_mydata", "tar_pheno", tar_covar_file = replaced_res, int_covar_index = 1, threads = 40)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
glm_df <- read.table(gweis_res$glm_file, header = TRUE, comment.char = "", stringsAsFactors = FALSE)

write.table(glm_df, file = "qq_gweis.phent.glm.linear", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

# Step 5: Munge for LDSC
cat("Step 5: Munging for LDSC...\n")
munge_res <- munge_ldsc_gcim(munge_path  = "/home/567/zw6700/ldsc/munge_sumstats.py", glm_file    = gweis_res, hm3_snplist = hm3_snps,
  python      = "/home/567/zw6700/anaconda3/envs/ldsc/bin/python")
# Step 6: LDSC heritability
cat("Step 6: Computing LDSC intercept...\n")
ldsc_res <- ldsc_h2_gcim(ldsc_path = "/home/567/zw6700/ldsc/ldsc.py", python = "/home/567/zw6700/anaconda3/envs/ldsc/bin/python", munged_sumstats = munge_res,
  ref_ld_chr = ld_scores)
# Step 7: Adjust Z-scores
cat("Step 7: Adjusting Z-scores...\n")
final_res <- gcim_z_adjust(glm_file = gweis_res, intercept_file = ldsc_res)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# fully adjusted outputs 
zdf <- read.table(final_res$output_file, header = TRUE, stringsAsFactors = FALSE)
file.copy(final_res$output_file, "gcim_gcim_z_qq_adjusted.txt", overwrite = TRUE)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
 ~~~

### 2. Quantitative Outcomes with Binary Exposures

~~~
#load library
library(GCIM.GWEIS.Z)
# Set plink path
plink_path <- "<plink_path>/plink2"
# Set munge_sumstats.py path
menge_path <- "<menge_path>/munge_sumstats.py"
# Set ldsc.py path
ldsc_path <- "<ldsc_path>/ldsc.py"
hm3_snps <- "<snplst_path>/w_hm3.snplist"
ld_scores <- "<ld_ref_path>/eur_w_ld_chr/"
# Step 1: Discovery GWAS
cat("Step 1: Running discovery GWAS...\n")
gwas_res <- b_gwas(plink_path = plink, "dis_mydata", "dis_cov", threads = 40)
# Step 2: Compute PRS
cat("Step 2: Computing PRS in target...\n")
prs_res <- prs_scores(plink_path = plink, "tar_mydata", score_file = gwas_res, threads = 40)

# Step 3: Replace exposure with PRS
cat("Step 3: Replacing exposure with PRS...\n")
replaced_res <- replace_covariate_with_prs("tar_covar", prs_file = prs_res, on_missing = "stop")

# Step 4: GWEIS
cat("Step 4: Running GWEIS...\n")
gweis_res <- q_gweis(plink_path = plink, "tar_mydata", "tar_pheno", tar_covar_file = replaced_res, int_covar_index = 1, threads = 40)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
glm_df <- read.table(gweis_res$glm_file, header = TRUE, comment.char = "", stringsAsFactors = FALSE)

write.table(glm_df, file = "qb_gweis.phent.glm.linear", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
# Step 5: Munge for LDSC
cat("Step 5: Munging for LDSC...\n")
munge_res <- munge_ldsc_gcim(munge_path  = "/home/567/zw6700/ldsc/munge_sumstats.py", glm_file    = gweis_res, hm3_snplist = hm3_snps, python      = "/home/567/zw6700/anaconda3/envs/ldsc/bin/python")
# Step 6: LDSC heritability
cat("Step 6: Computing LDSC intercept...\n")
ldsc_res <- ldsc_h2_gcim(ldsc_path = "/home/567/zw6700/ldsc/ldsc.py", python = "/home/567/zw6700/anaconda3/envs/ldsc/bin/python", munged_sumstats = munge_res, ref_ld_chr = ld_scores)

# Step 7: Adjust Z-scores
cat("Step 7: Adjusting Z-scores...\n")
final_res <- gcim_z_adjust(glm_file = gweis_res, intercept_file = ldsc_res)
#save the adjusted values 
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# fully adjusted outputs 
zdf <- read.table(final_res$output_file, header = TRUE, stringsAsFactors = FALSE)
file.copy(final_res$output_file, "gcim_gcim_z_qb_adjusted.txt", overwrite = TRUE)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
 ~~~

### 3. Binary Outcomes with Quantitative Exposures 

~~~
#load library
library(GCIM.GWEIS.Z)
# Set plink path
plink_path <- "<plink_path>/plink2"
# Set munge_sumstats.py path
menge_path <- "<menge_path>/munge_sumstats.py"
# Set ldsc.py path
ldsc_path <- "<ldsc_path>/ldsc.py"
hm3_snps <- "<snplst_path>/w_hm3.snplist"
ld_scores <- "<ld_ref_path>/eur_w_ld_chr/"

# 1) GWAS -> score file
# Step 1: Discovery GWAS
cat("Step 1: Running discovery GWAS...\n")
gwas_res <- q_gwas(plink_path = plink, "dis_mydata", "dis_pheno", threads = 40)

# Step 2: Compute PRS
cat("Step 2: Computing PRS in target...\n")
prs_res <- prs_scores(plink_path = plink, "tar_mydata", score_file = gwas_res, threads = 40)
# Step 3: Replace exposure with PRS
cat("Step 3: Replacing exposure with PRS...\n")
replaced_res <- replace_covariate_with_prs("tar_covar", prs_file = prs_res, on_missing = "stop")

# Step 4: GWEIS
cat("Step 4: Running GWEIS...\n")
gweis_res <- b_gweis(plink_path = plink, "tar_mydata", "tar_pheno", tar_covar_file = replaced_res, int_covar_index = 1, threads = 40)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
glm_df <- read.table(gweis_res$glm_file, header = TRUE, comment.char = "", stringsAsFactors = FALSE)

write.table(glm_df, file = "bq_gweis.phent.glm.linear", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
# Step 5: Munge for LDSC
cat("Step 5: Munging for LDSC...\n")
munge_res <- munge_ldsc_gcim(munge_path  = "/home/567/zw6700/ldsc/munge_sumstats.py", glm_file    = gweis_res, hm3_snplist = hm3_snps,
  python      = "/home/567/zw6700/anaconda3/envs/ldsc/bin/python")
# Step 6: LDSC heritability
cat("Step 6: Computing LDSC intercept...\n")
ldsc_res <- ldsc_h2_gcim(ldsc_path = "/home/567/zw6700/ldsc/ldsc.py", python = "/home/567/zw6700/anaconda3/envs/ldsc/bin/python",
  munged_sumstats = munge_res, ref_ld_chr = ld_scores)

# Step 7: Adjust Z-scores
cat("Step 7: Adjusting Z-scores...\n")
final_res <- gcim_z_adjust(glm_file = gweis_res, intercept_file = ldsc_res)
#save the adjusted values 
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# fully adjusted outputs 
zdf <- read.table(final_res$output_file, header = TRUE, stringsAsFactors = FALSE)
file.copy(final_res$output_file, "gcim_gcim_z_bq_adjusted.txt", overwrite = TRUE)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
 ~~~

### 4. Binary Outcomes with Binary Exposures

~~~
#load library
library(GCIM.GWEIS.Z)
# Set plink path
plink_path <- "<plink_path>/plink2"
# Set munge_sumstats.py path
menge_path <- "<menge_path>/munge_sumstats.py"
# Set ldsc.py path
ldsc_path <- "<ldsc_path>/ldsc.py"
hm3_snps <- "<snplst_path>/w_hm3.snplist"
ld_scores <- "<ld_ref_path>/eur_w_ld_chr/"

# Step 1: Discovery GWAS
cat("Step 1: Running discovery GWAS...\n")
gwas_res <- b_gwas(plink_path = plink, "dis_mydata", "dis_pheno", threads = 40)

# Step 2: Compute PRS
cat("Step 2: Computing PRS in target...\n")
prs_res <- prs_scores(plink_path = plink, "tar_mydata", score_file = gwas_res, threads = 40)

# Step 3: Replace exposure with PRS
cat("Step 3: Replacing exposure with PRS...\n")
replaced_res <- replace_covariate_with_prs("tar_covar", prs_file = prs_res, on_missing = "stop")

# Step 4: GWEIS
cat("Step 4: Running GWEIS...\n")
gweis_res <- b_gweis(plink_path = plink, "tar_mydata", "tar_pheno", tar_covar_file = replaced_res, int_covar_index = 1, threads = 40)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
glm_df <- read.table(gweis_res$glm_file, header = TRUE, comment.char = "", stringsAsFactors = FALSE)

write.table(glm_df, file = "bb_gweis.phent.glm.linear", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
# Step 5: Munge for LDSC
cat("Step 5: Munging for LDSC...\n")
munge_res <- munge_ldsc_gcim(munge_path  = "/home/567/zw6700/ldsc/munge_sumstats.py", glm_file    = gweis_res,
  hm3_snplist = hm3_snps, python      = "/home/567/zw6700/anaconda3/envs/ldsc/bin/python")
# Step 6: LDSC heritability
cat("Step 6: Computing LDSC intercept...\n")
ldsc_res <- ldsc_h2_gcim(ldsc_path = "/home/567/zw6700/ldsc/ldsc.py", python = "/home/567/zw6700/anaconda3/envs/ldsc/bin/python",
  munged_sumstats = munge_res, ref_ld_chr = ld_scores)

# Step 7: Adjust Z-scores
cat("Step 7: Adjusting Z-scores...\n")
final_res <- gcim_z_adjust(glm_file = gweis_res, intercept_file = ldsc_res)

#save the adjusted values 
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# fully adjusted outputs 
zdf <- read.table(final_res$output_file, header = TRUE, stringsAsFactors = FALSE)
file.copy(final_res$output_file, "gcim_gcim_z_bb_adjusted.txt", overwrite = TRUE)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
 ~~~

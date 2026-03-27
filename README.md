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
### Aditional softwares 
To run the GCIM-GWEIS-Z pipeline, the following external software must be installed and accessible from your system: 
   [Plink](https://www.cog-genomics.org/plink/2.0/) and 
   [LDSC](https://github.com/bulik/LDSC). 
## 2. Load the library


~~~
library(GCIM.GWEIS.Z)
~~~

## 3. Data Preparation 
 
To ensure robust and unbiased estimation of G×E interactions, the dataset must be divided into two independent, non-overlapping subsets: a discovery dataset and a target dataset. This split must be applied consistently across all inputs, including genotype data, phenotypes (outcomes), environmental exposures, and covariates.

The discovery dataset is used to estimate SNP effects, while the target dataset is used for downstream evaluation and GCIM-based inference.

Example datasets are provided in the data/ directory. These illustrate two causal directions:

***Body mass index (BMI) as the outcome and total bilirubin as the exposure***, and
***Total bilirubin as the outcome and BMI as the exposure***.

These examples are included to demonstrate how the input structure should be defined for both directions of analysis.

#### 3.1. Genotype data 
Genotype data must be provided in PLINK binary format, consisting of three files: `.bed`, `.bim`, and `.fam`.

Separate genotype datasets are required for the discovery and target samples. These datasets must contain non-overlapping individuals and be consistently aligned with the corresponding phenotype and covariate files.

Example files are provided in the data/ directory:

***Discovery dataset***
`mydata_dis.bed`
`mydata_dis.bim`
`mydata_dis.fam`
***Target dataset***
`mydata_tar.bed`
`mydata_tar.bim`
`mydata_tar.fam`
Ensure that standard quality control (QC) procedures (e.g., SNP missingness, minor allele frequency, and Hardy–Weinberg equilibrium filtering) are applied before analysis.
### 3.2. Input File Format (Discovery Dataset for GWAS/PRS Construction)
The discovery dataset is used to perform GWAS and construct polygenic risk scores (PRS) for the exposure variable under the proposed causal directions. 
Exposure + Covariates File
The input file must include the following columns: `tBil_dis_cov.txt`
`FID`
`IID`
`bilirubin`
`TDI, Age, Sex, ..., Covar_n`

If the exposure is binary, use PLINK’s default coding :

`1` = Control
`2` = Case

#### 3.3. Genome-wide environment interaction study (GWEIS)

The file should contain two separate `.txt` files for outcome ( `bmi_tar_out.txt`) and exposure with covariate files ('tBil_tar_cov.txt'): This is a .txt file containing the following columns in the specified order. Please note that the file should not have column headings. Therefore, the outcome file `bmi_tar_out.txt` will have the following essential column:

* FID

* IID

* BMI


'tBil_tar_cov.txt': This is a .txt file containing the following columns in the specified order. Note that the file should have no column heading. The exposure variable and the covariate that are used to adjust the data frame, as expressed in GxEprs. 

* FID

* IID
  
*  Bilirubin  

* TDI
* age
* sex
  
  .
  .
  .
  
* covar_n
  
***NB*** The reverse direction analysis follows the same input file format and structure as described above. The only difference is the switching of roles between the outcome and exposure variables. Specifically: The variable previously treated as the outcome is now used as the exposure, and the variable previously treated as the exposure is now used as the outcome.
Example
Proposed causal: Exposure = Total bilirubin and Outcome = BMI
Reverse direction: Exposure = BMI and Outcome = Total bilirubin

All other components (genotype data, covariates, file structure, and coding conventions) remain unchanged.

## 4. Analysis pipeline
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
# Step 4: GWEIS (This is a step 1 for conventional GWEIS)
cat("Step 4: Running GWEIS...\n")
gweis_res <- q_gweis(plink_path = plink, "tar_mydata", "tar_pheno", tar_covar_file = replaced_res, int_covar_index = 1, threads = 40)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
      ###save the GCIM-GWEIS
#glm_df <- read.table(gweis_res$glm_file, header = TRUE, comment.char = "", stringsAsFactors = FALSE)

#write.table(glm_df, file = "qq_gweis.phent.glm.linear", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
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
# fully adjusted GCIM_GWEIS-Z outputs 
zdf <- read.table(final_res$output_file, header = TRUE, stringsAsFactors = FALSE)
file.copy(final_res$output_file, "gcim_gcim_z_qq_adjusted.txt", overwrite = TRUE)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
 ~~~
###Result 


<img width="1900" height="1139" alt="image" src="https://github.com/user-attachments/assets/2d30e953-e85c-4804-a654-16144045ed82" />


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




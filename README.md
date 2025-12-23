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
plink_path <- "<plink_path>/plink2"
ldsc_path <- "<ldsc_path>/ldsc.py"
munge_path <- "<munge_path>/munge_sumstats.py"
plink_path <- "<plink_path>/plink2"
~~~
library(GCIM-GWEIS-Z)
~~~
# Quantitative trait

~~~
#load library
library(GCIM-GWEIS-Z)
# specify the file paths

# compute GWAS of the Exposure variable
gwas <- q_gwas(plink_path, dis_snp, qp_dis_cov)
#compute PRS of Exposure variables
prs <- prs_scores(plink_path, tar_snp, score_file)
# Merged PRS values with covariate files
replace_covariate_with_prs(
  qp_dis_cov,
  prs_file,
  out_name = "qp_dis_cov_prs.txt"
)
#compute GWEIS analyses
gweis <- q_gweis(
  plink_path,
  dis_snp,
  qp_dis_phen,
  qp_dis_cov,
  int_covar_index = 1,
  out_file = "gcim_prs_int.txt"
)
# compute heritability for the interaction components using LDSC
ldsc <- run_ldsc_gcim(
  ldsc_path,
  munge_path,
  gcim_file,
  hm3_snplist,
  ref_ld_chr,
  chunksize = 1e+05
)
# Adjusting z-score values of GWEIS results.
gcim-gweis-z <-gcim_z_adjust(
  glm_file,
  intercept_file,
  int_term = "ADDxCOVAR1",
  out_file = "gcim_z_adjusted.txt"
)
write.table(gcim-gweis-z, "GCIM-GWEIS-Z.txt", quote = F, row.names = F, col.names = T, sep = " ")
 

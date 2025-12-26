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
library(GCIM-GWEIS-Z)
~~~

## 3. Analysis pipeline
### 1. Quantitative Outcomes with Quantitative Exposures 

~~~
#load library
library(GCIM-GWEIS-Z)
# Set plink path
plink_path <- "<plink_path>/plink2"
# Set munge_sumstats.py path
menge_path <- "<menge_path>/munge_sumstats.py"
# Set ldsc.py path
ldsc_path <- "<ldsc_path>/ldsc.py"

# 1) GWAS -> score file
g <- q_gwas(plink_path, dis_mydata, dis_cov_file)

# 2) PRS in target
p <- prs_scores(plink_path, tar_mydata, g$score_file)

# 3) Make covariate file where COVAR1 = PRS (column 3)
cov_prs <- replace_covariate_with_prs(tar_covar_file, p$prs_file)

# 4) GWEIS (quant or binary)
gw <- q_gweis(plink_path, tar_mydata, tar_pheno_file, cov_prs$output_file, int_covar_index = 1)

# 5) Munge for LDSC
m <- munge_ldsc_gcim(munge_path, gw, gw$interaction_term, hm3_snplist)

# 6) LDSC h2 / intercept
h2 <- ldsc_h2_gcim(ldsc_path, m, ref_ld_chr, w_ld_chr = ref_ld_chr)

# 7) Z adjust
zad <- gcim_z_adjust(gw, h2, int_term = gw$interaction_term)
#save the adjusted values 
write.table(zad, "GCIM-GWEIS-Z_qq.txt", quote = F, row.names = F, col.names = T, sep = " ")
 ~~~
### 2. Quantitative Outcomes with Binary Exposures

~~~
#load library
library(GCIM-GWEIS-Z)

# Set plink path
plink_path <- "<plink_path>/plink2"
# Set munge_sumstats.py path
menge_path <- "<menge_path>/munge_sumstats.py"
# Set ldsc.py path
ldsc_path <- "<ldsc_path>/ldsc.py"

# 1) GWAS -> score file
g <- b_gwas(plink_path, dis_mydata, dis_cov_file)

# 2) PRS in target
p <- prs_scores(plink_path, tar_mydata, g$score_file)

# 3) Make covariate file where COVAR1 = PRS (column 3)
cov_prs <- replace_covariate_with_prs(tar_covar_file, p$prs_file)

# 4) GWEIS (quant or binary)
gw <- q_gweis(plink_path, tar_mydata, tar_pheno_file, cov_prs$output_file, int_covar_index = 1)

# 5) Munge for LDSC
m <- munge_ldsc_gcim(munge_path, gw, gw$interaction_term, hm3_snplist)

# 6) LDSC h2 / intercept
h2 <- ldsc_h2_gcim(ldsc_path, m, ref_ld_chr, w_ld_chr = ref_ld_chr)

# 7) Z adjust
zad <- gcim_z_adjust(gw, h2, int_term = gw$interaction_term)
#save the adjusted values 
write.table(zad, "GCIM-GWEIS-Z_qb.txt", quote = F, row.names = F, col.names = T, sep = " ")
 ~~~

### 3. Binary Outcomes with Quantitative Exposures 

~~~
#load library
library(GCIM-GWEIS-Z)

# Set plink path
plink_path <- "<plink_path>/plink2"
# Set munge_sumstats.py path
menge_path <- "<menge_path>/munge_sumstats.py"
# Set ldsc.py path
ldsc_path <- "<ldsc_path>/ldsc.py"

# 1) GWAS -> score file
g <- q_gwas(plink_path, dis_mydata, dis_cov_file)

# 2) PRS in target
p <- prs_scores(plink_path, tar_mydata, g$score_file)

# 3) Make covariate file where COVAR1 = PRS (column 3)
cov_prs <- replace_covariate_with_prs(tar_covar_file, p$prs_file)

# 4) GWEIS (quant or binary)
gw <- b_gweis(plink_path, tar_mydata, tar_pheno_file, cov_prs$output_file, int_covar_index = 1)

# 5) Munge for LDSC
m <- munge_ldsc_gcim(munge_path, gw, gw$interaction_term, hm3_snplist)

# 6) LDSC h2 / intercept
h2 <- ldsc_h2_gcim(ldsc_path, m, ref_ld_chr, w_ld_chr = ref_ld_chr)

# 7) Z adjust
zad <- gcim_z_adjust(gw, h2, int_term = gw$interaction_term)
#save the adjusted values 
write.table(zad, "GCIM-GWEIS-Z_bq.txt", quote = F, row.names = F, col.names = T, sep = " ")
 ~~~

### 4. Binary Outcomes with Binary Exposures

~~~
#load library
library(GCIM-GWEIS-Z)

# Set plink path
plink_path <- "<plink_path>/plink2"
# Set munge_sumstats.py path
menge_path <- "<menge_path>/munge_sumstats.py"
# Set ldsc.py path
ldsc_path <- "<ldsc_path>/ldsc.py"

# 1) GWAS -> score file
g <- b_gwas(plink_path, dis_mydata, dis_cov_file)

# 2) PRS in target
p <- prs_scores(plink_path, tar_mydata, g$score_file)

# 3) Make covariate file where COVAR1 = PRS (column 3)
cov_prs <- replace_covariate_with_prs(tar_covar_file, p$prs_file)

# 4) GWEIS (quant or binary)
gw <- b_gweis(plink_path, tar_mydata, tar_pheno_file, cov_prs$output_file, int_covar_index = 1)

# 5) Munge for LDSC
m <- munge_ldsc_gcim(munge_path, gw, gw$interaction_term, hm3_snplist)

# 6) LDSC h2 / intercept
h2 <- ldsc_h2_gcim(ldsc_path, m, ref_ld_chr, w_ld_chr = ref_ld_chr)

# 7) Z adjust
zad <- gcim_z_adjust(gw, h2, int_term = gw$interaction_term)
#save the adjusted values 
write.table(zad, "GCIM-GWEIS-Z_bb.txt", quote = F, row.names = F, col.names = T, sep = " ")
 ~~~

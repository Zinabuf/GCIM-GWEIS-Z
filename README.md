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
#### Additional Resources for LD Score Regression
**HapMap3 SNP list**:
`w_hm3.snplist` (used for SNP filtering and harmonization)  and 
**European LD scores**:
`eur_w_ld_chr/` (used as the LD reference panel)

## 2. Load the library

~~~
library(GCIM.GWEIS.Z)
~~~

## 3. Data Preparation 
 
The dataset must be divided into two independent, non-overlapping subsets: a discovery dataset and a target dataset. This split must be applied consistently across all inputs, including genotype data, phenotypes (outcomes), environmental exposures, and covariates.

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
gwas_res <- q_gwas(plink_path = plink, "mydata.dis", dis_cov_file = "tBil_dis_cov.txt", threads = 40)
# Step 2: Compute PRS
cat("Step 2: Computing PRS in target...\n")
prs_res <- prs_scores(plink_path = plink, "mydata.tar", score_file = gwas_res, threads = 40)

# Step 3: Replace exposure with PRS
cat("Step 3: Replacing exposure with PRS...\n")
replaced_res <- replace_covariate_with_prs("tBil_tar_cov.txt", prs_file = prs_res, on_missing = "stop")
# Step 4: GWEIS (This is a step 1 for conventional GWEIS)
cat("Step 4: Running GWEIS...\n")
gweis_res <- q_gweis(plink_path = plink, "mydata.tar", "bmi_tar_out.txt", tar_covar_file = replaced_res, int_covar_index = 1, threads = 40)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
      ###save the GCIM-GWEIS
#glm_df <- read.table(gweis_res$glm_file, header = TRUE, comment.char = "", stringsAsFactors = FALSE)

#write.table(glm_df, file = "tBil_bmi_gcim-gweis.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

# Step 5: Munge for LDSC
cat("Step 5: Munging for LDSC...\n")
munge_res <- munge_ldsc_gcim(munge_path  = "/<path>/ldsc/munge_sumstats.py", glm_file    = gweis_res, hm3_snplist = hm3_snps,
  python      = "/<path>/anaconda3/envs/ldsc/bin/python")
# Step 6: LDSC heritability
cat("Step 6: Computing LDSC intercept...\n")
ldsc_res <- ldsc_h2_gcim(ldsc_path = "/<path>/ldsc/ldsc.py", python = "/<path>/anaconda3/envs/ldsc/bin/python", munged_sumstats = munge_res,
  ref_ld_chr = ld_scores)
# Step 7: Adjust Z-scores
cat("Step 7: Adjusting Z-scores...\n")
final_res <- gcim_z_adjust(glm_file = gweis_res, intercept_file = ldsc_res)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# fully adjusted GCIM_GWEIS-Z outputs 
zdf <- read.table(final_res$output_file, header = TRUE, stringsAsFactors = FALSE)
file.copy(final_res$output_file, "tBil_bmi_gcim-gweis-z.txt", overwrite = TRUE)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
 ~~~
### Result 
GCIM-GWEIS
~~~
X.CHROM	POS	ID	REF	ALT	PROVISIONAL_REF.	A1	OMITTED	A1_FREQ	TEST	OBS_CT	BETA	SE	T_STAT	P	ERRCODE
1	5348270	rs4417082	G	A	Y	A	G	0.42228	ADD	193	0.0673094	0.101462	0.663392	0.507958	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	PRS	193	-0.117611	0.128938	-0.912155	0.362949	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	tdi	193	0.00814524	0.0759784	0.107205	0.91475	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	age	193	-0.0464425	0.0756185	-0.614169	0.539906	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	sex	193	0.0479832	0.072578	0.661126	0.509406	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	PC1	193	-0.0517115	0.0776805	-0.665695	0.506488	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	PC2	193	-0.0809538	0.0740406	-1.09337	0.275743	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	PC3	193	-0.13183	0.0738943	-1.78403	0.076161	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	PC4	193	-0.189733	0.102401	-1.85285	0.0655981	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	PC5	193	0.112447	0.104436	1.07671	0.283101	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	PC6	193	0.0333456	0.0766078	0.435277	0.663901	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	PC7	193	-0.0193206	0.0817745	-0.236267	0.813503	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	PC8	193	-0.0164211	0.0831785	-0.197419	0.84373	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	PC9	193	0.0409719	0.082082	0.499159	0.618298	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	PC10	193	0.1279	0.0799992	1.59877	0.111686	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	center	193	0.0554347	0.0747765	0.741338	0.459488	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	batch	193	0.0299986	0.0728423	0.41183	0.680971	.
1	5348270	rs4417082	G	A	Y	A	G	0.42228	ADDxPRS	193	0.0832006	0.108047	0.770043	0.442318	.
~~~

GCIM-GWEUIS-Z

~~~
#CHROM	POS	ID	REF	ALT	PROVISIONAL_REF?	A1	OMITTED	A1_FREQ	TEST	OBS_CT	BETA	SE	T_STAT	P	ERRCODE	z_original	z_adj_int	p_original	p_value_int_adj	ldsc_intercept_used
1	5348270	rs4417082	G	A	Y	A	G	0.42228	ADDxPRS	193	0.0832006	0.108047	0.770043	0.442318	.	0.770043	0.793858762886598	0.442318	0.427277658606437	0.9409
1	8287391	rs11121137	G	A	Y	A	G	0.234848	ADDxPRS	198	-0.147242	0.109512	-1.34452	0.180482	.	-1.34452	-1.38610309278351	0.180482	0.165715426128863	0.9409
1	8407287	rs10399665	T	C	Y	C	T	0.4575	ADDxPRS	200	0.0710143	0.10361	0.685402	0.493967	.	0.685402	0.7066	0.493967	0.47981508894707	0.9409
1	12642243	rs12029485	T	G	Y	G	T	0.0351759	ADDxPRS	199	-0.28364	0.392133	-0.723326	0.470418	.	-0.723326	-0.745696907216495	0.470418	0.455850527893872	0.9409
1	16612193	rs12739564	T	G	Y	G	T	0.010101	ADDxPRS	198	-1.63974	1.24591	-1.31609	0.189826	.	-1.31609	-1.35679381443299	0.189826	0.174846733475931	0.9409
1	17550601	rs3003482	A	G	Y	G	A	0.1625	ADDxPRS	200	0.0168795	0.150742	0.111976	0.910966	.	0.111976	0.115439175257732	0.910966	0.908097029138002	0.9409
1	22653550	rs4433361	C	T	Y	T	C	0.347716	ADDxPRS	197	-0.055233	0.117466	-0.470203	0.638785	.	-0.470203	-0.484745360824742	0.638785	0.627856985392955	0.9409
1	30198580	rs2376723	A	G	Y	G	A	0.268229	ADDxPRS	192	-0.174419	0.107968	-1.61548	0.10803	.	-1.61548	-1.66544329896907	0.10803	0.0958243470239232	0.9409
1	30274404	rs16832918	G	A	Y	A	G	0.223077	ADDxPRS	195	-0.0455193	0.136091	-0.334477	0.738418	.	-0.334477	-0.344821649484536	0.738418	0.730228456238774	0.9409
1	30661242	rs4245643	A	G	Y	G	A	0.43	ADDxPRS	200	-0.0431465	0.121007	-0.356563	0.721834	.	-0.356563	-0.367590721649485	0.721834	0.713178434225341	0.9409
1	37041069	rs4653197	C	T	Y	T	C	0.2	ADDxPRS	195	-0.0495329	0.139871	-0.354132	0.723664	.	-0.354132	-0.365084536082474	0.723664	0.715048303822747	0.9409
1	37621304	rs4652937	C	T	Y	T	C	0.333333	ADDxPRS	195	0.0408107	0.107049	0.381234	0.703489	.	0.381234	0.393024742268041	0.703489	0.694301207147922	0.9409
1	39952822	rs1569053	C	T	Y	T	C	0.1125	ADDxPRS	200	0.0204852	0.179559	0.114086	0.909296	.	0.114086	0.117614432989691	0.909296	0.906373168453006
~~~

Number of Significant SNPs

~~~
Z-score adjustment completed successfully!
  Correction factor (sqrt(intercept)): 0.970000
  Valid SNPs: 1000

Significance counts (p < 0.05):
  Original:  59 / 1000 (5.90%)
  Adjusted:  67 / 1000 (6.70%)

Bonferroni-corrected (p < 5.00e-05):
  Original:  0
  Adjusted:  0
~~~


### 2. Binary Outcomes with Quantitative Exposures 

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
munge_res <- munge_ldsc_gcim(munge_path  = "/<path>/ldsc/munge_sumstats.py", glm_file    = gweis_res, hm3_snplist = hm3_snps,
  python      = "/<path>/anaconda3/envs/ldsc/bin/python")
# Step 6: LDSC heritability
cat("Step 6: Computing LDSC intercept...\n")
ldsc_res <- ldsc_h2_gcim(ldsc_path = "/<path>/ldsc/ldsc.py", python = "/<path>/anaconda3/envs/ldsc/bin/python",
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
### Result
<img width="1929" height="1073" alt="image" src="https://github.com/user-attachments/assets/66ee2335-0f5f-41cd-bf75-1ba6a7a5546d" />




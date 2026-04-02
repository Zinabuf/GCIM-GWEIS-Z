---
# <h1 align="center">GCIM-GWEIS-Z</h1>
---

The Genetic Causality Inference Model (GCIM) for genome-wide-by-environment interaction studies (GWEIS), incorporating adjusted z-scores (Z), is a statistical framework designed to infer the causal direction of SNP-by-environment interactions and to correct inflation in GxE effects arising from heteroscedasticity.
## Authors: Zinabu Fentaw and S.Hong Lee

   
# 1. Package installation 
From GitHub 

~~~
library(devtools)
install_github("Zinabuf/GCIM-GWEIS-Z")
~~~
### Aditional softwares 
To run the GCIM-GWEIS-Z pipeline, the following external software must be installed and accessible from your system: 
   [Plink](https://www.cog-genomics.org/plink/2.0/) and 
   [LDSC](https://github.com/bulik/LDSC). 
### Additional Resources for LD Score Regression
**HapMap3 SNP list**:
`w_hm3.snplist` (used for SNP filtering and harmonization)  and 
**European LD scores**:
`eur_w_ld_chr/` (used as the LD reference panel)

# 2. Load the library

~~~
library(GCIM.GWEIS.Z)
~~~

# 3. Data Preparation 
 
The dataset must be divided into two independent, non-overlapping subsets: a PRS training sample and an analysis sample. This split must be applied consistently across all inputs, including genotype data, phenotypes (outcomes), environmental exposures, and covariates.

The PRS training sample was used to obtain SNP effect size estimates from GWAS summary statistics for polygenic risk score construction, whereas the independent analysis sample was used to perform GWEIS and evaluate genotype-by-environment(GXE) interactions, including inference on causal direction and assessment of test statistic inflation.

Example datasets are provided in the `data/` directory. All datasets—including genotype data, phenotypes (outcome and exposure), and covariates are simulated. These examples illustrate the causal directional analysis of G×E interaction effects:

 **Trait 1 as the outcome and Trait 2 as the exposure**, and
 **Trait 2 as the outcome and Trait 1 as the exposure**.

They are intended to demonstrate how to correctly specify the input structure for each causal direction. In each case, the exposure variable should be defined within a dedicated data frame, along with any covariates included for adjustment in the model.
 
## 3.1. Genotype data 
Genotype data must be provided in PLINK binary format, consisting of three files: `.bed`, `.bim`, and `.fam`.

Separate genotype datasets are required for the PRS training and analysis samples. These datasets must contain non-overlapping individuals and be consistently aligned with the corresponding phenotype and covariate files.

Example files are provided in the data/ directory:

***PRS training sample***
`mydata.PRStrained.bed`
`mydata.PRStrained.bim`
`mydata.PRStrained.fan`

***Analysis sample***
`mydata.analysis.bed`
`mydata.analysis.bim`
`mydata.analysis.fam`
## 3.2. Input File Format (Discovery Dataset for GWAS/PRS Construction)
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

## 3.3. Genome-wide environment interaction study (GWEIS)

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

# 4. Analysis pipeline
## 1. Quantitative Outcomes with Quantitative Exposures 

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
## Result 
Here are the first few lines of a typical **GCIM-GWEIS** output file for quantitative outcome `trait1_out_q_gcim_gweis.txt`

~~~
X.CHROM POS     ID      REF     ALT     PROVISIONAL_REF.        A1      OMITTED A1_FREQ TEST    OBS_CT  BETA    SE      T_STAT  P       ERRCODE
1       1000    rs1     A       G       Y       G       A       0.181875        ADD     800     0.148359        0.0848825       1.74782 0.0808831       .
1       1000    rs1     A       G       Y       G       A       0.181875        PRS     800     -0.0778738      0.0557167       -1.39767        0.162603        .
1       1000    rs1     A       G       Y       G       A       0.181875        TDI     800     -0.0122875      0.0463571       -0.265062       0.791031        .
1       1000    rs1     A       G       Y       G       A       0.181875        age     800     -0.0675298      0.0461932       -1.4619 0.144166        .
1       1000    rs1     A       G       Y       G       A       0.181875        PC1     800     0.00015977      0.0462788       0.00345234      0.997246        .
1       1000    rs1     A       G       Y       G       A       0.181875        PC2     800     -0.0251955      0.046216        -0.545168       0.585792        .
1       1000    rs1     A       G       Y       G       A       0.181875        PC3     800     -0.0136874      0.0461746       -0.296427       0.766982        .
1       1000    rs1     A       G       Y       G       A       0.181875        ADDxPRS 800     0.195433        0.0864984       2.25938 0.0241313       .
1       2000    rs23    A       G       Y       G       A       0.34    ADD     800     -0.0301563      0.0684743       -0.440403       0.659766        .
1       2000    rs23    A       G       Y       G       A       0.34    PRS     800     -0.00848817     0.0665136       -0.127616       0.898486        .
1       2000    rs23    A       G       Y       G       A       0.34    TDI     800     -0.0185266      0.0464972       -0.398446       0.690409        .
1       2000    rs23    A       G       Y       G       A       0.34    age     800     -0.0662399      0.0464724       -1.42536        0.154448        .
1       2000    rs23    A       G       Y       G       A       0.34    PC1     800     -0.0010281      0.0465056       -0.022107       0.982368        .
1       2000    rs23    A       G       Y       G       A       0.34    PC2     800     -0.0302483      0.0465483       -0.649826       0.515993        .
1       2000    rs23    A       G       Y       G       A       0.34    PC3     800     -0.0189549      0.0463756       -0.408727       0.682851        .
1       2000    rs23    A       G       Y       G       A       0.34    ADDxPRS 800     -0.00462024     0.0685233       -0.0674257      0.94626 .
~~~

Here are the first few lines of a typical **GCIM-GWEIS-Z** output file for quantitative outcome `trait1_out_q_gcim_gweis.txt`

~~~
#CHROM	POS	ID	REF	ALT	PROVISIONAL_REF?	A1	OMITTED	A1_FREQ	TEST	OBS_CT	BETA	SE	T_STAT	P	ERRCODE	z_original	z_adj_int	p_original	p_value_int_adj	ldsc_intercept_used
1	1000	rs1	A	G	Y	G	A	0.181875	ADDxPRS	800	0.195433	0.0864984	2.25938	0.0241313	.	2.25938	2.57647605235805	0.0241313	0.00998131192882362	0.769
1	2000	rs23	A	G	Y	G	A	0.34	ADDxPRS	800	-0.00462024	0.0685233	-0.0674257	0.94626	.	-0.0674257	-0.0768886603242828	0.94626	0.938712118591888	0.769
1	3000	rs45	A	G	Y	G	A	0.125	ADDxPRS	800	-0.120493	0.0680116	-1.77165	0.0768379	.	-1.77165	-2.02029485883744	0.0768379	0.0433528117744768	0.769
1	4000	rs67	A	G	Y	G	A	0.41	ADDxPRS	800	-0.107577	0.0679464	-1.58326	0.113763	.	-1.58326	-1.80546498360453	0.113763	0.0710019502377452	0.769
1	5000	rs89	A	G	Y	G	A	0.441875	ADDxPRS	800	-0.0254022	0.0642314	-0.39548	0.692596	.	-0.39548	-0.450984229827015	0.692596	0.652000915047852	0.769
1	6000	rs111	A	G	Y	G	A	0.49	ADDxPRS	800	0.0283261	0.0463104	0.611659	0.540939	.	0.611659	0.697503193667852	0.540939	0.48548794143712	0.769
1	7000	rs133	A	G	Y	G	A	0.3275	ADDxPRS	800	-0.0380365	0.0490007	-0.776246	0.437836	.	-0.776246	-0.885189401401591	0.437836	0.376054487633956	0.769
1	8000	rs155	A	G	Y	G	A	0.115	ADDxPRS	800	-0.118236	0.096108	-1.23024	0.218972	.	-1.23024	-1.40289986573882	0.218972	0.160646700930315	0.769
1	9000	rs177	A	G	Y	G	A	0.213125	ADDxPRS	800	0.186653	0.0843303	2.21335	0.0271581	.	2.21335	2.52398590342779	0.0271581	0.0116032581361642	0.769
1	10000	rs199	A	G	Y	G	A	0.10875	ADDxPRS	800	-0.0187808	0.080353	-0.233729	0.815256	.	-0.233729	-0.266532044738643	0.815256	0.789829483523949	0.769
1	11000	rs221	A	G	Y	G	A	0.1725	ADDxPRS	800	-0.0987807	0.0574429	-1.71963	0.0858903	.	-1.71963	-1.96097403443266	0.0858903	0.0498820519342425	0.769
1	12000	rs243	A	G	Y	G	A	0.11875	ADDxPRS	800	-0.0974253	0.0758968	-1.28366	0.199638	.	-1.28366	-1.46381717522946	0.199638	0.143243908344006	0.769
1	13000	rs265	A	G	Y	G	A	0.275625	ADDxPRS	800	0.0438517	0.0696596	0.629514	0.529194	.	0.629514	0.717864080245078	0.529194	0.472841097468128	0.769
1	14000	rs287	A	G	Y	G	A	0.0775	ADDxPRS	800	-0.111224	0.0933534	-1.19143	0.233841	.	-1.19143	-1.35864301846566	0.233841	0.174259736486416	0.769
1	15000	rs309	A	G	Y	G	A	0.26375	ADDxPRS	800	-0.0462584	0.0509923	-0.907165	0.364596	.	-0.907165	-1.03448242351326	0.364596	0.30091068951762	0.769
1	16000	rs331	A	G	Y	G	A	0.3175	ADDxPRS	800	-0.0809935	0.0691342	-1.17154	0.241735	.	-1.17154	-1.33596152678148	0.241735	0.181561848782453	0.769
~~~

Number of Significant SNPs

~~~
Z-score adjustment completed successfully!
  Correction factor (sqrt(intercept)): 0.970000
  Valid SNPs: 1000

Significance counts (p < 0.05):
  Original(GCIM-GWEIS):  59 / 1000 (5.90%)
  Adjusted(GCIM-GWEIS-Z):  67 / 1000 (6.70%)

Bonferroni-corrected (p < 5.00e-05):
  Original(GCIM-GWEIS):  0
  Adjusted(GCIM-GWEIS-Z):  0
~~~


## 2. Binary Outcomes with Quantitative Exposures 

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
gwas_res <- q_gwas(plink_path = plink, "mydata.dis", "tBil_dis_cov.txt", threads = 40)

# Step 2: Compute PRS
cat("Step 2: Computing PRS in target...\n")
prs_res <- prs_scores(plink_path = plink, "mydata.tar", score_file = gwas_res, threads = 40)
# Step 3: Replace exposure with PRS
cat("Step 3: Replacing exposure with PRS...\n")
replaced_res <- replace_covariate_with_prs("tBil_tar_cov.txt", prs_file = prs_res, on_missing = "stop")

# Step 4: GWEIS
cat("Step 4: Running GWEIS...\n")
gweis_res <- b_gweis(plink_path = plink, "mydata.tar", "bmi_tar_out_b.txt", tar_covar_file = replaced_res, int_covar_index = 1, threads = 40)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
glm_df <- read.table(gweis_res$glm_file, header = TRUE, comment.char = "", stringsAsFactors = FALSE)

write.table(glm_df, file = "tBil_bmi_gcim-gweis.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
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
file.copy(final_res$output_file, "tBil_bmi_gcim-gweis-z.txt", overwrite = TRUE)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
 ~~~
## Result

Here are the first few lines of a typical **GCIM-GWEIS** output file `trait1_out_b_gcim_gweis.txt`
 
 ~~~
X.CHROM POS     ID      REF     ALT     PROVISIONAL_REF.        A1      OMITTED A1_FREQ FIRTH.  TEST    OBS_CT  OR      LOG.OR._SE      Z_STAT  P       ERRCODE
1       1000    rs1     A       G       Y       G       A       0.181875        N       ADD     800     1.15489 0.140329        1.02621 0.30479 .
1       1000    rs1     A       G       Y       G       A       0.181875        N       PRS     800     0.9199  0.0946198       -0.882375       0.377574        .
1       1000    rs1     A       G       Y       G       A       0.181875        N       TDI     800     1.0035  0.0778832       0.0448172       0.964253        .
1       1000    rs1     A       G       Y       G       A       0.181875        N       age     800     0.923445        0.0778339       -1.02326        0.306186        .
1       1000    rs1     A       G       Y       G       A       0.181875        N       PC1     800     0.977338        0.0778073       -0.294605       0.768296        .
1       1000    rs1     A       G       Y       G       A       0.181875        N       PC2     800     0.987061        0.0775306       -0.16798        0.866599        .
1       1000    rs1     A       G       Y       G       A       0.181875        N       PC3     800     1.03138 0.0774197       0.399113        0.68981 .
1       1000    rs1     A       G       Y       G       A       0.181875        N       ADDxPRS 800     1.18503 0.143839        1.18024 0.237905        .
1       2000    rs23    A       G       Y       G       A       0.34    N       ADD     800     0.934197        0.114914        -0.592338       0.553624        .
1       2000    rs23    A       G       Y       G       A       0.34    N       PRS     800     1.02179 0.11006 0.195822        0.84475 .
1       2000    rs23    A       G       Y       G       A       0.34    N       TDI     800     0.995121        0.0775005       -0.0631132      0.949676        .
1       2000    rs23    A       G       Y       G       A       0.34    N       age     800     0.923729        0.0777545       -1.02034        0.307565        .
1       2000    rs23    A       G       Y       G       A       0.34    N       PC1     800     0.977041        0.0777822       -0.298609       0.765238        .
1       2000    rs23    A       G       Y       G       A       0.34    N       PC2     800     0.985724        0.0777548       -0.184927       0.853286        .
1       2000    rs23    A       G       Y       G       A       0.34    N       PC3     800     1.02575 0.0774383       0.328319        0.742671        .
1       2000    rs23    A       G       Y       G       A       0.34    N       ADDxPRS 800     0.935743        0.114551        -0.579779       0.562064        .
~~~

Here are the first few lines of a typical **GCIM-GWEIS-Z** output file `trait1_out_b_gcim_gweis-z.txt`

~~~
#CHROM	POS	ID	REF	ALT	PROVISIONAL_REF?	A1	OMITTED	A1_FREQ	FIRTH?	TEST	OBS_CT	OR	LOG(OR)_SE	Z_STAT	P	ERRCODE	z_original	z_adj_int	p_original	p_value_int_adj	ldsc_intercept_used
1	1000	rs1	A	G	Y	G	A	0.181875	N	ADDxPRS	800	1.18503	0.143839	1.18024	0.237905	.	1.18024	1.36537484636713	0.237905	0.172135265534288	0.7472
1	2000	rs23	A	G	Y	G	A	0.34	N	ADDxPRS	800	0.935743	0.114551	-0.579779	0.562064	.	-0.579779	-0.67072431289559	0.562064	0.502396171410051	0.7472
1	3000	rs45	A	G	Y	G	A	0.125	N	ADDxPRS	800	0.862357	0.138155	-1.07189	0.28377	.	-1.07189	-1.24002884504208	0.28377	0.214964725292673	0.7472
1	4000	rs67	A	G	Y	G	A	0.41	N	ADDxPRS	800	0.815581	0.113836	-1.79078	0.0733293	.	-1.79078	-2.07168539227389	0.0733293	0.0382947905081261	0.7472
1	5000	rs89	A	G	Y	G	A	0.441875	N	ADDxPRS	800	0.968699	0.107914	-0.294689	0.768232	.	-0.294689	-0.340914515777371	0.768232	0.733167936732633	0.7472
1	6000	rs111	A	G	Y	G	A	0.49	N	ADDxPRS	800	0.938498	0.0775154	-0.818859	0.412867	.	-0.818859	-0.947306887854458	0.412867	0.343482422035134	0.7472
1	7000	rs133	A	G	Y	G	A	0.3275	N	ADDxPRS	800	0.889865	0.0811446	-1.43799	0.150437	.	-1.43799	-1.66355603549064	0.150437	0.0962011839371125	0.7472
1	8000	rs155	A	G	Y	G	A	0.115	N	ADDxPRS	800	0.859721	0.167999	-0.899689	0.368286	.	-0.899689	-1.04081604601878	0.368286	0.297960931448276	0.7472
1	9000	rs177	A	G	Y	G	A	0.213125	N	ADDxPRS	800	1.45284	0.152704	2.44602	0.0144443	.	2.44602	2.82970767107617	0.0144443	0.00465905523120323	0.7472
1	10000	rs199	A	G	Y	G	A	0.10875	N	ADDxPRS	800	1.21976	0.133462	1.48847	0.136627	.	1.48847	1.72195443093955	0.136627	0.0850777736247182	0.7472
1	11000	rs221	A	G	Y	G	A	0.1725	N	ADDxPRS	800	0.907518	0.0968739	-1.00174	0.316471	.	-1.00174	-1.15887497339508	0.316471	0.246507151474062	0.7472
1	12000	rs243	A	G	Y	G	A	0.11875	N	ADDxPRS	800	0.924044	0.125084	-0.631539	0.527688	.	-0.631539	-0.730603491747318	0.527688	0.465021378629281	0.7472
1	13000	rs265	A	G	Y	G	A	0.275625	N	ADDxPRS	800	1.03162	0.11522	0.270163	0.787035	.	0.270163	0.312541317544808	0.787035	0.754629167662689	0.7472
1	14000	rs287	A	G	Y	G	A	0.0775	N	ADDxPRS	800	0.967925	0.147627	-0.220832	0.825223	.	-0.220832	-0.255472156572347	0.825223	0.798358446105053	0.7472
1	15000	rs309	A	G	Y	G	A	0.26375	N	ADDxPRS	800	1.06068	0.0847563	0.695097	0.486994	.	0.695097	0.804131328869769	0.486994	0.421321132990375	0.7472
1	16000	rs331	A	G	Y	G	A	0.3175	N	ADDxPRS	800	0.822538	0.119603	-1.63341	0.102383	.	-1.63341	-1.88963001406878	0.102383	0.058807460110183	0.7472
1	17000	rs353	A	G	Y	G	A	0.276875	N	ADDxPRS	800	0.988598	0.12406	-0.0924312	0.926355	.	-0.0924312	-0.106930145986859	0.926355	0.914844397255351	0.7472
1	18000	rs375	A	G	Y	G	A	0.298125	N	ADDxPRS	800	0.996369	0.112922	-0.0322169	0.974299	.	-0.0322169	-0.0372705084456766	0.974299	0.970269319999561	0.7472
1	19000	rs397	A	G	Y	G	A	0.16625	N	ADDxPRS	800	0.994236	0.0923127	-0.0626177	0.950071	.	-0.0626177	-0.0724400397523922	0.950071	0.942251721362135	0.7472
~~~
 
Number of Significant SNPs

~~~
Z-score adjustment completed successfully!
  Correction factor (sqrt(intercept)): 0.864407
  Valid SNPs: 1000

Significance counts (p < 0.05):
  Original:  54 / 1000 (5.40%)
  Adjusted:  88 / 1000 (8.80%)

Bonferroni-corrected for number of SNPs(p < 5.00e-05):
  Original:  0
  Adjusted:  0

Genome-wide significance threshold(p < 5e-8):
  Original_genomewide_bonf:  0
  Adjusted_genomewide_bonf:  0
~~~

All results presented above were generated using **BMI as the outcome** and **total bilirubin as the exposure**. Example input files and code for other outcome-exposure combinations are also provided in the `data/` directory. For comparison of your output with the toy example results, please refer to the corresponding `.Rout` files.


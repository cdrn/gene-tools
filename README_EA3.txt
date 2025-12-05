*******************************************************************************
*** 1. File Contents
*******************************************************************************

*** Column Names: 

MarkerName: SNP rs number.
CHR: chromosome number.
POS: base pair position.
A1: effect allele.
A2: other allele.
EAF: A1 frequency in 1000 Genomes Phase 3 sample (CEU, GBR and TSI individuals).
Beta_*: Standardized regression coefficient, i.e. per-allele effect size on the phenotype that has been standardized to have unit variance.
SE_*: standard error of Beta
Pval_*: Nominal p-value of the null hypothesis that the coefficient is equal to zero.

*** SNP Selection: 

To limit the possibilities of identifiability, significant digits for betas and standard errors are limited to five decimal points and no sample allele frequencies are provided (the reported allele frequencies are calculated using data on European-ancestry individuals who contributed to the 1000 Genomes Project - Phase 3). Association results are only provided for SNPs that pass standard quality-control filters described in the SI of Lee et al. (2018). We have also imposed a sample-size filter of 0.5*n_max, where n_max is the maximum sample size in each analysis.

*******************************************************************************
*** 2. All-SNP Summary Statistics
*******************************************************************************

For all SNPs that survive the QC and sample-size filters, we provide association results from the following GWAS meta-analyses:

1. GWAS_EA_excl23andMe.txt - Educational attainment (EA) GWAS meta-analysis of all discovery cohorts except 23andMe. Sample Size = 766,345.

2. GWAS_CP_all.txt - Cognitive performance (CP) GWAS meta-analysis of all discovery cohorts. Sample Size = 257,828.

*******************************************************************************
*** 3. 10K-SNP Summary Statistics
*******************************************************************************

Files 3-11 below are restricted to SNPs satisfying at least one of the following conditions:

    (i) Lead SNPs (p-value < 5e-8) in the analyses of Cognitive Performance (CP), self-reported math ability (MA) or highest math class completed (HM).
    (ii) Lead SNPs for at least one phenotype in the MTAG analysis of EA, CP, HM, and MA.
    (iii) SNPs that remain after applying clumping algorithm to summary statistics from the primary EA GWAS (all cohorts) at a p-value threshold of 1e-4.

Conditions (i) to (iii) were selected to ensure the final number of SNP-level summary statistics from analyses that include 23andMe does not exceed 10,000. Below, we refer to this subset of SNPs as our "10K SNPs". Users are advised that some of the 10K SNPs are in linkage equilibrium.

3. GWAS_EA.to10k.txt - Educational attainment meta-analysis of all discovery cohorts for 10K SNPs. Sample Size = 1,131,881.

4. GWAS_CP.to10K.txt - Cognitive performance meta-analysis of all discovery cohorts for 10K SNPs. Sample Size = 257,828.

5. GWAS_HM.to10K.txt - Highest-level math class completed GWAS for 10K SNPs. Sample Size = 430,439.

6. GWAS_MA.to10K.txt - Self-reported math ability GWAS for 10K SNPs. Sample Size = 564,692.

7. MTAG_EA.to10K.txt - EA results from MTAG analysis with CP, HM, and MA. The MTAG analysis was restricted to the 10K SNPs.

8. MTAG_CP.to10K.txt - CP results from MTAG analysis with EA, HM, and MA. The MTAG analysis was restricted to the 10K SNPs.

9. MTAG_HM.to10K.txt - HM results from MTAG analysis with CP, EA, and MA. The MTAG analysis was restricted to the 10K SNPs.

10. MTAG_MA.to10K.txt - MA results from MTAG analysis with CP, HM, and MA. The MTAG analysis was restricted to the 10K SNPs.

11. COMBINED.to10K.txt - Combined results from files 3-10, with the corresponding columns of each result suffixed by analysis type and trait (e.g., "Beta_GWAS_HM").

12. CAVIARBF.to10K.txt - Results from CAVIARBF analysis for 10K SNPs. PIP column gives posterior probability SNP is causal. NA means that the SNP was not included in this analysis. [Added on Sep 07, 2018]

For additional details, please see the Supplementary Note accompanying Lee et al. (2018).

Reference

Lee et al. (2018). Gene discovery and polygenic prediction from a 1.1-million-person GWAS of educational attainment. Nature Genetics, 50 (8), 1112-1121. doi: 10.1038/s41588-018-0147-3
# Set working directory
setwd("../UKB/IS-EWAS&MR/MR_PRESSO")
workdir <- getwd()


# Load packages
pkgs <- c(
    "dplyr", "tidyr", "readxl", "readr", "data.table", "tidyverse", "survival",
    "gwasglue", "gwasvcf", "MRInstruments", "VariantAnnotation",
    "TwoSampleMR", "MRPRESSO", "RadialMR", "mr.raps", "MendelianRandomization"
)
inst <- lapply(pkgs, require, character.only = TRUE)


# Load data
plink_freq <- read.table(
    "../UKB/IS-EWAS&MR/data_archived/Plink/your_output.frq",
    header = TRUE
)


# Check the status of the API
ieugwasr::api_status()
ieugwasr::get_opengwas_jwt()
ieugwasr::user()


# 1 IEU: GWAS ID----

## 1.1 Exposure data P<5*10-8 ----
IEU_ID <- read.csv("../UKB/IS-EWAS&MR/data_archived/GWAS_Exposure_ID.csv") %>%
    filter(GWAS_ID != "") %>%
    dplyr::pull(GWAS_ID) %>%
    as.character()


# Create empty data frame
result_combine_f <- data.frame()


# MR-PRESSO
for (i in IEU_ID) {
    print(paste("Current loop index:", i))


    # Load Exposure Data
    exposure_dat <- read.csv(paste0(
        "../UKB/IS-EWAS&MR/MR_10-8/Exposure_Resources/", i, "_clumped.csv"
    ))
    trait1_exposure_data_clumped <- exposure_dat

    if (!is.null(exposure_dat)) {
        # Load Outcome Data
        outcome_dat <- TwoSampleMR::extract_outcome_data(
            snps = trait1_exposure_data_clumped$SNP,
            outcomes = "ebi-a-GCST006908",
            proxies = TRUE, # LD proxies
            # By default if a particular requested SNP is not present in the outcome GWAS then a SNP (proxy) that is in LD with the requested SNP (target) will be searched for instead.
            # LD proxies are defined using 1000 genomes European sample data.
            # The effect of the proxy SNP on the outcome is returned, along with the proxy SNP, the effect allele of the proxy SNP, and the corresponding allele (in phase) for the target SNP.
            rsq = 0.8,
            align_alleles = 1, # Try to align tag alleles to target alleles
            palindromes = 1, # If TRUE, palindromic SNPs will be flipped to match the effect allele of the target SNP
            maf_threshold = 0.3, # Minimum MAF for proxy SNPs
            opengwas_jwt = ieugwasr::get_opengwas_jwt(),
            splitsize = 10000, # Number of SNPs to process at a time
            proxy_splitsize = 500 # Number of proxy SNPs to process at a time
        )

        if (!is.null(outcome_dat)) {
            # Harmonise data
            trait1_trait2_dat <- TwoSampleMR::harmonise_data(
                exposure_dat = trait1_exposure_data_clumped,
                outcome_dat = outcome_dat,
                action = 2
            )


            # Perform MR
            trait1_trait2_results <- TwoSampleMR::mr(
                trait1_trait2_dat,
                parameters = TwoSampleMR::default_parameters()
            )
            trait1_trait2_results_with_or <- TwoSampleMR::generate_odds_ratios(
                trait1_trait2_results
            )
            result_combine <- filter(
                trait1_trait2_results_with_or,
                method == "Inverse variance weighted"
            )


            # MR-PRESSO needs >3 SNPs
            if (nrow(trait1_trait2_dat) > 2) {
                trait1_trait2_run <- TwoSampleMR::run_mr_presso(
                    dat = trait1_trait2_dat,
                    NbDistribution = 1000, # 8000
                    SignifThreshold = 0.05
                )

                result_combine$MR_PRESSO_Global_pval <-
                    trait1_trait2_run[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
                result_combine$Main_MR_results.Effect_Raw <-
                    trait1_trait2_run[[1]]$`Main MR results`$`Causal Estimate`[1]
                result_combine$Main_MR_results.Sd_Raw <-
                    trait1_trait2_run[[1]]$`Main MR results`$Sd[1]
                result_combine$Main_MR_results.pval_Raw <-
                    trait1_trait2_run[[1]]$`Main MR results`$`P-value`[1]
                result_combine$Main_MR_results.Effect_Corre <-
                    trait1_trait2_run[[1]]$`Main MR results`$`Causal Estimate`[2]
                result_combine$Main_MR_results.Sd_Corre <-
                    trait1_trait2_run[[1]]$`Main MR results`$Sd[2]
                result_combine$Main_MR_results.pval_Corre <-
                    trait1_trait2_run[[1]]$`Main MR results`$`P-value`[2]
                result_combine$MR_PRESSO.Results.Distortion.pval <-
                    trait1_trait2_run[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue
                if (is.null(result_combine$MR_PRESSO.Results.Distortion.pval)) {
                    result_combine$MR_PRESSO.Results.Distortion.pval <- NA
                }
                result_combine$MR_PRESSO_Global_pval <-
                    as.character(result_combine$MR_PRESSO_Global_pval)
                result_combine$MR_PRESSO.Results.Distortion.pval <-
                    as.character(result_combine$MR_PRESSO.Results.Distortion.pval)

                result_combine_f <- dplyr::bind_rows(
                    result_combine_f, result_combine
                )
                write.csv(
                    result_combine_f,
                    "../UKB/IS-EWAS&MR/MR_PRESSO/Results_10-8/result_combine_run.csv",
                    row.names = FALSE
                )
            }
        }
    }
}


# 2 FinnGen: PhenoCode----

## 2.1 Exposure data P<5*10-8----
PhenoCode <- read.csv(
    "../UKB/IS-EWAS&MR/data_archived/FinnGen/PhenoCode_5e-8.csv"
) %>%
    filter(PhenoCode != "") %>%
    dplyr::pull(PhenoCode) %>%
    as.character()

# Create empty data frame
result_combine_f <- data.frame()

# MR-PRESSO
for (i in PhenoCode) {
    print(paste("Current loop index:", i))

    # Load Exposure Data
    exposure_dat <- read.csv(paste0(
        "../UKB/IS-EWAS&MR/MR/EWAS 68-Forward-FinnGen/MR_10-8/Exposure_Resources/",
        i,
        "_clumped.csv"
    ))
    trait1_exposure_data_clumped <- exposure_dat

    if (!is.null(exposure_dat)) {
        # Load Outcome Data
        outcome_dat <- TwoSampleMR::extract_outcome_data(
            snps = trait1_exposure_data_clumped$SNP,
            outcomes = "ebi-a-GCST006908",
            proxies = TRUE, # LD proxies
            # By default if a particular requested SNP is not present in the outcome GWAS then a SNP (proxy) that is in LD with the requested SNP (target) will be searched for instead.
            # LD proxies are defined using 1000 genomes European sample data.
            # The effect of the proxy SNP on the outcome is returned, along with the proxy SNP, the effect allele of the proxy SNP, and the corresponding allele (in phase) for the target SNP.
            rsq = 0.8,
            align_alleles = 1, # Try to align tag alleles to target alleles
            palindromes = 1, # If TRUE, palindromic SNPs will be flipped to match the effect allele of the target SNP
            maf_threshold = 0.3, # Minimum MAF for proxy SNPs
            opengwas_jwt = ieugwasr::get_opengwas_jwt(),
            splitsize = 10000, # Number of SNPs to process at a time
            proxy_splitsize = 500 # Number of proxy SNPs to process at a time
        )

        if (!is.null(outcome_dat)) {
            # Harmonise data
            trait1_trait2_dat <- TwoSampleMR::harmonise_data(
                exposure_dat = trait1.exposure_data_clumped,
                outcome_dat = outcome_dat,
                action = 2
            )


            # Perform MR
            trait1_trait2_results <- TwoSampleMR::mr(
                trait1_trait2_dat,
                parameters = TwoSampleMR::default_parameters()
            )
            trait1_trait2_results_with_or <-
                TwoSampleMR::generate_odds_ratios(trait1_trait2_results)
            result_combine <- filter(
                trait1_trait2_results_with_or,
                method == "Inverse variance weighted"
            )


            # MR-PRESSO needs >3 SNPs
            if (nrow(trait1_trait2_dat) > 2) {
                run_mr_presso <- TwoSampleMR::run_mr_presso(
                    dat = trait1_trait2_dat,
                    NbDistribution = 8000, # 8000
                    SignifThreshold = 0.05
                )

                result_combine$MR_PRESSO_Global_pval <-
                    run_mr_presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
                result_combine$Main_MR_results.Effect_Raw <-
                    run_mr_presso[[1]]$`Main MR results`$`Causal Estimate`[1]
                result_combine$Main_MR_results.Sd_Raw <-
                    run_mr_presso[[1]]$`Main MR results`$Sd[1]
                result_combine$Main_MR_results.pval_Raw <-
                    run_mr_presso[[1]]$`Main MR results`$`P-value`[1]
                result_combine$Main_MR_results.Effect_Corre <-
                    run_mr_presso[[1]]$`Main MR results`$`Causal Estimate`[2]
                result_combine$Main_MR_results.Sd_Corre <-
                    run_mr_presso[[1]]$`Main MR results`$Sd[2]
                result_combine$Main_MR_results.pval_Corre <-
                    run_mr_presso[[1]]$`Main MR results`$`P-value`[2]
                result_combine$MR_PRESSO.Results.Distortion.pval <-
                    run_mr_presso[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue
                if (is.null(result_combine$MR_PRESSO.Results.Distortion.pval)) {
                    result_combine$MR_PRESSO.Results.Distortion.pval <- NA
                }
                result_combine$MR_PRESSO_Global_pval <-
                    as.character(result_combine$MR_PRESSO_Global_pval)
                result_combine$MR_PRESSO.Results.Distortion.pval <-
                    as.character(result_combine$MR_PRESSO.Results.Distortion.pval)

                result_combine_f <- dplyr::bind_rows(
                    result_combine_f, result_combine
                )
                write.csv(
                    result_combine_f,
                    "../UKB/IS-EWAS&MR/MR_PRESSO/Results_10-8/result_combine_run.csv",
                    row.names = FALSE
                )
            }
        }
    }
}

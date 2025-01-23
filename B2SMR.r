# Set working directory
setwd("../UKB/IS-EWAS&MR/MR")
workdir <- getwd()


# Load packages
pkgs <- c(
    "dplyr", "tidyr", "readxl", "tidyverse", "data.table", "survival",
    "gwasglue", "gwasvcf", "MRInstruments", "VariantAnnotation",
    "TwoSampleMR", "simex", "RadialMR", "mr.raps", "MendelianRandomization"
)
inst <- lapply(pkgs, require, character.only = TRUE)


# Create empty data frame
result_combine_f <- data.frame()
df_norun <- data.frame()


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
IEU_ID <- read.csv("../UKB/IS-EWAS&MR/data_archived/GWAS_Exposure_ID.csv") %>%
    filter(GWAS_ID != "") %>%
    dplyr::pull(GWAS_ID) %>%
    as.character()

for (i in IEU_ID) {
    print(paste("Current loop index:", i))

    # Load Exposure Data
    exposure_dat <- TwoSampleMR::extract_instruments(
        outcomes = i,
        p1 = 5e-08, # P-value threshold for keeping a SNP
        clump = FALSE, # Whether or not to return independent SNPs only
        p2 = 5e-08, # P-value threshold for clumping
        r2 = 0.001, # The maximum LD R² allowed between returned SNPs
        kb = 10000, # The distance in which to search for LD R² values
        opengwas_jwt = ieugwasr::get_opengwas_jwt(),
        force_server = TRUE
    )

    if (!is.null(exposure_dat)) {
        trait1_exposure_data <- exposure_dat
        write.csv(
            trait1_exposure_data,
            paste0("../UKB/IS-EWAS&MR/MR/Exposure_Resources/", i, ".csv"),
            row.names = FALSE
        )


        # Clumping by running local LD operations
        trait1_exposure_data <- dplyr::rename(
            trait1_exposure_data,
            rsid = SNP,
            pval = pval.exposure,
            trait_id = id.exposure
        )
        trait1_exposure_data_clumped <- ieugwasr::ld_clump(
            dat = trait1_exposure_data,
            clump_kb = 10000,
            clump_r2 = 0.001,
            clump_p = 0.99,
            pop = "EUR",
            bfile = "../UKB/IS-EWAS&MR/plink_1kg.v3/EUR",
            plink_bin = genetics.binaRies::get_plink_binary()
        )

        trait1_exposure_data_clumped <- dplyr::rename(
            trait1_exposure_data_clumped,
            SNP = rsid,
            pval.exposure = pval,
            id.exposure = id
        )
        write.csv(
            trait1.exposure_data_clumped,
            paste0(
                "../UKB/IS-EWAS&MR/MR/Exposure_Resources/",
                i,
                "_clumped.csv"
            ),
            row.names = FALSE
        )


        # Load Outcome Data
        outcome_dat <- TwoSampleMR::extract_outcome_data(
            snps = trait1.exposure_data_clumped$SNP,
            outcomes = "ebi-a-GCST006908",
            proxies = TRUE, # LD proxies
            # By default if a particular requested SNP is not present in the outcome GWAS then a SNP (proxy) that is in LD with the requested SNP (target) will be searched for instead.
            # LD proxies are defined using 1000 genomes European sample data.
            # The effect of the proxy SNP on the outcome is returned, along with the proxy SNP, the effect allele of the proxy SNP, and the corresponding allele (in phase) for the target SNP.
            rsq = 0.8, # numeric value of minimum rsq to find a proxy. Default is 0.8, minimum is 0.6
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
            write.csv(
                trait1_trait2_dat,
                paste0(
                    "../UKB/IS-EWAS&MR/MR/Harmonise_Data/",
                    i,
                    "_ebi-a-GCST006908_Harmonise.csv"
                ),
                row.names = FALSE
            )


            # Perform MR
            trait1_trait2_results <- TwoSampleMR::mr(
                trait1_trait2_dat,
                parameters = TwoSampleMR::default_parameters()
            )
            trait1_trait2_results_with_or <-
                TwoSampleMR::generate_odds_ratios(trait1_trait2_results)


            # Depends on SNP counts
            if (nrow(trait1_trait2_dat) > 1) {
                # Turn the results into one row
                result_combine <- dplyr::filter(
                    trait1_trait2_results_with_or,
                    method == "Inverse variance weighted"
                )
                result_MREGGER <- dplyr::filter(
                    trait1_trait2_results_with_or,
                    method == "MR Egger"
                )
                result_WMe <- dplyr::filter(
                    trait1_trait2_results_with_or,
                    method == "Weighted median"
                )
                result_WMo <- dplyr::filter(
                    trait1_trait2_results_with_or,
                    method == "Weighted mode"
                )
                result_SM <- dplyr::filter(
                    trait1_trait2_results_with_or,
                    method == "Simple mode"
                )
                result_combine <- dplyr::left_join(
                    result_combine, result_MREGGER,
                    by = c("id.exposure", "id.outcome", "outcome", "exposure")
                )
                result_combine <- dplyr::left_join(
                    result_combine, result_WMe,
                    by = c("id.exposure", "id.outcome", "outcome", "exposure")
                )
                result_combine <- dplyr::left_join(
                    result_combine, result_WMo,
                    by = c("id.exposure", "id.outcome", "outcome", "exposure")
                )
                result_combine <- dplyr::left_join(
                    result_combine, result_SM,
                    by = c("id.exposure", "id.outcome", "outcome", "exposure")
                )


                # MR-EGGER INTERCEPT
                pleiotropy_test_result <-
                    TwoSampleMR::mr_pleiotropy_test(trait1_trait2_dat)
                result_combine$egger_intercept <-
                    pleiotropy_test_result$egger_intercept
                result_combine$egger_intercept.pval <-
                    pleiotropy_test_result$pval

                # MR-EGGER/MR-EGGER heterogeneity need at least three SNPs
                if (nrow(trait1_trait2_dat) > 2) {
                    # MR-EGGER-SIMEX
                    BetaYG <- trait1_trait2_dat$beta.outcome
                    BetaXG <- trait1_trait2_dat$beta.exposure
                    seBetaYG <- trait1_trait2_dat$se.outcome
                    seBetaXG <- trait1_trait2_dat$se.exposure

                    Fit2 <- lm(
                        BetaYG ~ BetaXG,
                        weights = 1 / seBetaYG^2,
                        x = TRUE,
                        y = TRUE
                    )
                    mod_sim <- simex::simex(
                        Fit2,
                        SIMEXvariable = "BetaXG",
                        measurement.error = seBetaXG,
                        B = 1000,
                        fitting.method = "quadratic",
                        asymptotic = FALSE
                    )
                    summary(mod_sim)
                    simex_beta <- summary(mod_sim)[[1]]$jackknife[2]
                    simex_se <- summary(mod_sim)[[1]]$jackknife[4]
                    simex_p <- summary(mod_sim)[[1]]$jackknife[8]

                    result_combine$isq <- TwoSampleMR::Isq(
                        trait1_trait2_results$b,
                        trait1_trait2_results$se
                    )
                    result_combine$simex_beta <- simex_beta
                    result_combine$simex_se <- simex_se
                    result_combine$simex_p <- simex_p


                    # Heterogeneity statistics
                    mr_heterogeneity_data <-
                        TwoSampleMR::mr_heterogeneity(
                            trait1_trait2_dat,
                            parameters = TwoSampleMR::default_parameters(),
                            method_list = subset(
                                TwoSampleMR::mr_method_list(),
                                heterogeneity_test & use_by_default
                            )$obj
                        )

                    result_combine$MR_Egger.Q <-
                        mr_heterogeneity_data$Q[
                            mr_heterogeneity_data$method == "MR Egger"
                        ]
                    result_combine$MR_Egger.Q_df <-
                        mr_heterogeneity_data$Q_df[
                            mr_heterogeneity_data$method == "MR Egger"
                        ]
                    result_combine$MR_Egger.Q_pval <-
                        mr_heterogeneity_data$Q_pval[
                            mr_heterogeneity_data$method == "MR Egger"
                        ]

                    result_combine$Inverse_variance_weighted.Q <-
                        mr_heterogeneity_data$Q[
                            mr_heterogeneity_data$method ==
                                "Inverse variance weighted"
                        ]
                    result_combine$Inverse_variance_weighted.Q_df <-
                        mr_heterogeneity_data$Q_df[
                            mr_heterogeneity_data$method ==
                                "Inverse variance weighted"
                        ]
                    result_combine$Inverse_variance_weighted.Q_pval <-
                        mr_heterogeneity_data$Q_pval[
                            mr_heterogeneity_data$method ==
                                "Inverse variance weighted"
                        ]
                } else {
                    # MR-EGGER-SIMEX
                    result_combine$isq <- NA
                    result_combine$isq <-
                        as.numeric(result_combine$isq)
                    result_combine$simex.beta <- NA
                    result_combine$simex.beta <-
                        as.numeric(result_combine$simex.beta)
                    result_combine$simex.se <- NA
                    result_combine$simex.se <-
                        as.numeric(result_combine$simex.se)
                    result_combine$simex.p <- NA
                    result_combine$simex.p <-
                        as.numeric(result_combine$simex.p)


                    # Heterogeneity statistics
                    mr_heterogeneity <- TwoSampleMR::mr_heterogeneity(
                        trait1_trait2_dat,
                        parameters = TwoSampleMR::default_parameters(),
                        method_list = subset(
                            TwoSampleMR::mr_method_list(),
                            heterogeneity_test & use_by_default
                        )$obj
                    )

                    result_combine$MR_Egger.Q <- NA
                    result_combine$MR_Egger.Q <-
                        as.numeric(result_combine$MR_Egger.Q)
                    result_combine$MR_Egger.Q_df <- NA
                    result_combine$MR_Egger.Q_df <-
                        as.numeric(result_combine$MR_Egger.Q_df)
                    result_combine$MR_Egger.Q_pval <- NA
                    result_combine$MR_Egger.Q_pval <-
                        as.numeric(result_combine$MR_Egger.Q_pval)

                    result_combine$Inverse_variance_weighted.Q <-
                        mr_heterogeneity$Q[
                            mr_heterogeneity$method ==
                                "Inverse variance weighted"
                        ]
                    result_combine$Inverse_variance_weighted.Q_df <-
                        mr_heterogeneity$Q_df[
                            mr_heterogeneity$method ==
                                "Inverse variance weighted"
                        ]
                    result_combine$Inverse_variance_weighted.Q_pval <-
                        mr_heterogeneity$Q_pval[
                            mr_heterogeneity$method ==
                                "Inverse variance weighted"
                        ]
                }


                # Generate reports
                newpath <- paste0(
                    "../UKB/IS-EWAS&MR/MR/Reports/", i, "_ebi-a-GCST006908"
                )
                dir.create(newpath)
                TwoSampleMR::mr_report(
                    dat = trait1_trait2_dat,
                    output_path = newpath,
                    output_type = "html",
                    author = "Manqi Zheng",
                    study = "Two Sample MR"
                )
            } else {
                # Turn the results into one row
                result_combine <- dplyr::filter(
                    trait1_trait2_results_with_or,
                    method == "Wald ratio"
                )
                result_MREGGER <- dplyr::filter(
                    trait1_trait2_results_with_or,
                    method == "MR Egger"
                )
                result_WMe <- dplyr::filter(
                    trait1_trait2_results_with_or,
                    method == "Weighted median"
                )
                result_WMo <- dplyr::filter(
                    trait1_trait2_results_with_or,
                    method == "Weighted mode"
                )
                result_SM <- dplyr::filter(
                    trait1_trait2_results_with_or,
                    method == "Simple mode"
                )
                result_combine <- dplyr::left_join(
                    result_combine, result_MREGGER,
                    by = c("id.exposure", "id.outcome", "outcome", "exposure")
                )
                result_combine <- dplyr::left_join(
                    result_combine, result_WMe,
                    by = c("id.exposure", "id.outcome", "outcome", "exposure")
                )
                result_combine <- dplyr::left_join(
                    result_combine, result_WMo,
                    by = c("id.exposure", "id.outcome", "outcome", "exposure")
                )
                result_combine <- dplyr::left_join(
                    result_combine, result_SM,
                    by = c("id.exposure", "id.outcome", "outcome", "exposure")
                )

                # MR-EGGER INTERCEPT
                result_combine$egger_intercept <- NA
                result_combine$egger_intercept <-
                    as.numeric(result_combine$egger_intercept)
                result_combine$egger_intercept.pval <- NA
                result_combine$egger_intercept.pval <-
                    as.numeric(result_combine$egger_intercept.pval)


                # MR-EGGER-SIMEX
                result_combine$isq <- NA
                result_combine$isq <- as.numeric(result_combine$isq)
                result_combine$simex.beta <- NA
                result_combine$simex.beta <-
                    as.numeric(result_combine$simex.beta)
                result_combine$simex.se <- NA
                result_combine$simex.se <- as.numeric(result_combine$simex.se)
                result_combine$simex.p <- NA
                result_combine$simex.p <- as.numeric(result_combine$simex.p)


                # Heterogeneity statistics
                result_combine$MR_Egger.Q <- NA
                result_combine$MR_Egger.Q <-
                    as.numeric(result_combine$MR_Egger.Q)
                result_combine$MR_Egger.Q_df <- NA
                result_combine$MR_Egger.Q_df <-
                    as.numeric(result_combine$MR_Egger.Q_df)
                result_combine$MR_Egger.Q_pval <- NA
                result_combine$MR_Egger.Q_pval <-
                    as.numeric(result_combine$MR_Egger.Q_pval)

                result_combine$Inverse_variance_weighted.Q <- NA
                result_combine$Inverse_variance_weighted.Q <-
                    as.numeric(result_combine$Inverse_variance_weighted.Q)
                result_combine$Inverse_variance_weighted.Q_df <- NA
                result_combine$Inverse_variance_weighted.Q_df <-
                    as.numeric(result_combine$Inverse_variance_weighted.Q_df)
                result_combine$Inverse_variance_weighted.Q_pval <- NA
                result_combine$Inverse_variance_weighted.Q_pval <-
                    as.numeric(result_combine$Inverse_variance_weighted.Q_pval)


                # Generate path but no report
                newpath <- paste0(
                    "../UKB/IS-EWAS&MR/MR/Reports/", i, "_ebi-a-GCST006908/"
                )
                dir.create(newpath)
            }


            # Other statistics
            lor <- trait1_trait2_dat$beta.exposure
            af <- dplyr::filter(
                plink_freq,
                SNP %in% trait1_trait2_dat$SNP
            )$MAF # SNP allele frequencey derived from 1kg plink bfile

            trait1_trait2_mr_steiger <- TwoSampleMR::mr_steiger2(
                r_exp = TwoSampleMR::get_r_from_pn(
                    trait1_trait2_dat$pval.exposure,
                    trait1_trait2_dat$samplesize.exposure
                ),
                r_out = TwoSampleMR::get_r_from_lor(
                    lor,
                    af,
                    ncase,
                    ncontrol,
                    ISprevalence,
                    model = "logit",
                    correction = FALSE
                ),
                n_exp = trait1_trait2_dat$samplesize.exposure,
                n_out = trait1_trait2_dat$samplesize.outcome
            )

            result_combine$steiger_test_p <-
                trait1_trait2_mr_steiger$steiger_test
            result_combine$causal_dir <-
                trait1_trait2_mr_steiger$correct_causal_direction
            result_combine$steiger_test_p_adj <-
                trait1_trait2_mr_steiger$steiger_test_adj
            result_combine$exp_r2 <-
                trait1_trait2_mr_steiger$r2_exp
            exp_r2 <- trait1_trait2_mr_steiger$r2_exp
            result_combine$f_stat <- (
                trait1_trait2_dat$samplesize.exposure[1] -
                    dim(trait1_trait2_dat)[1] - 1
            ) / (dim(trait1_trait2_dat)[1]) * exp_r2 / (1 - exp_r2)


            # Calculate power
            OR <- trait1_trait2_results_with_or$or[
                trait1_trait2_results_with_or$method ==
                    "Inverse variance weighted" |
                    trait1_trait2_results_with_or$method == "Wald ratio"
            ]
            rsq <- exp.R2
            result_combine$power <- pnorm(
                sqrt(
                    n * rsq * (ratio / (1 + ratio)) * (1 / (1 + ratio))
                ) * OR - qnorm(1 - sig / 2)
            )

            write.csv(
                result_combine,
                paste0(
                    "../UKB/IS-EWAS&MR/MR/Results_Combine_Individual/",
                    i,
                    "_ebi-a-GCST006908", ".csv"
                ),
                row.names = FALSE
            )


            # Single SNP analysis
            setwd(newpath)
            single_snp_analysis <- TwoSampleMR::mr_singlesnp(
                trait1_trait2_dat,
                parameters = TwoSampleMR::default_parameters(),
                single_method = "mr_wald_ratio",
                all_method = c("mr_ivw", "mr_egger_regression")
            )

            single_snp_analysis_with_or <-
                TwoSampleMR::generate_odds_ratios(
                    single_snp_analysis
                )
            write.csv(
                single_snp_analysis_with_or,
                paste0("SSA_", i, "_ebi-a-GCST006908", ".csv"),
                row.names = FALSE
            )


            # Leave-one-out analysis
            mr_leaveoneout <- TwoSampleMR::mr_leaveoneout(
                trait1_trait2_dat,
                parameters = TwoSampleMR::default_parameters()
            )
            mr_leaveoneout_with_or <-
                TwoSampleMR::generate_odds_ratios(
                    mr_leaveoneout
                )

            write.csv(
                mr_leaveoneout_with_or,
                paste0("LOO_", i, "_ebi-a-GCST006908", ".csv"),
                row.names = FALSE
            )
            write.csv(
                result_combine_f,
                "../UKB/IS-EWAS&MR/MR/Results_Combine_Individual/Result_combine_run.csv",
                row.names = FALSE
            )


            # Combine results
            result_combine_f <- dplyr::bind_rows(
                result_combine_f, result_combine
            )

            write.csv(
                result_combine_f,
                "../UKB/IS-EWAS&MR/MR/Results_Combine_Individual/Result_combine_run.csv",
                row.names = FALSE
            )
        } else {
            # Record the condition where no SNP was found in the outcome dataset
            reason_norun <- "nosigsnp_in_outcome"
            df_norun_add <- data.frame(
                fieldid = i,
                reason_norun = reason_norun
            )
            df_norun <- rbind.data.frame(df_norun, df_norun_add)

            write.csv(
                df_norun,
                "../UKB/IS-EWAS&MR/MR/Results_Combine_Individual/NOsigSNP_ls.csv",
                row.names = FALSE
            )
        }
    } else {
        # Record the condition where no SNP was found in the exposure dataset
        reason_norun <- "NOsigSNP_in_Exposure"
        df_norun_add <- data.frame(
            fieldid = i,
            reason_norun = reason_norun
        )
        df_norun <- rbind.data.frame(df_norun, df_norun_add)

        write.csv(
            df_norun,
            "../UKB/IS-EWAS&MR/MR/Results_Combine_Individual/NOsigSNP_ls.csv",
            row.names = FALSE
        )
    }
}


# 2 FinnGen: PhenoCode----

## 2.1 Generate the exposure data----

# Common data processing
data <- data.table::fread(exposureFile)
output <- subset(data, pval < 5e-08)
output$ncase <- ncase
output$ncontrol <- ncontrol
output$samplesize <- output$ncase + output$ncontrol
nobs <- nrow(output)
if (nobs > 2) {
    output_dir <- "../UKB/IS-EWAS&MR/data_archived/FinnGen/finngen_R12_5e-8"
    output_file <- file.path(output_dir, paste0(exposureName, ".csv"))
    write.csv(output, file = output_file, row.names = FALSE)
} else if (nobs == 1 || nobs == 2) {
    output_dir_5e8 <- "../UKB/IS-EWAS&MR/data_archived/FinnGen/finngen_R12_5e-8"
    output_file_5e8 <- file.path(output_dir_5e8, paste0(exposureName, ".csv"))
    write.csv(output, file = output_file_5e8, row.names = FALSE)

    output_5e6 <- subset(data, pval < 5e-06)
    output_5e6$ncase <- ncase
    output_5e6$ncontrol <- ncontrol
    output_5e6$samplesize <- output_5e6$ncase + output_5e6$ncontrol
    nobs_5e6 <- nrow(output_5e6)
    if (nobs_5e6 > 0) {
        output_dir_5e6 <-
            "../UKB/IS-EWAS&MR/data_archived/FinnGen/finngen_R12_5e-6"
        output_file_5e6 <- file.path(
            output_dir_5e6, paste0(exposureName, ".csv")
        )
        write.csv(output_5e6, file = output_file_5e6, row.names = FALSE)
    } else if (nobs_5e6 == 0) {
        message("No significant observations found. No 5e-6 data saved.")
    }
} else if (nobs == 0) {
    message("No significant observations found. No 5e-8 data saved.")

    output_5e6 <- subset(data, pval < 5e-06)
    output_5e6$ncase <- ncase
    output_5e6$ncontrol <- ncontrol
    output_5e6$samplesize <- output_5e6$ncase + output_5e6$ncontrol
    nobs_5e6 <- nrow(output_5e6)
    if (nobs_5e6 > 0) {
        output_dir_5e6 <-
            "../UKB/IS-EWAS&MR/data_archived/FinnGen/finngen_R12_5e-6"
        output_file_5e6 <- file.path(
            output_dir_5e6, paste0(exposureName, ".csv")
        )
        write.csv(output_5e6, file = output_file_5e6, row.names = FALSE)
    } else if (nobs_5e6 == 0) {
        message("No significant observations found. No 5e-6 data saved.")
    }
}


## 2.2 Exposure data P<5*10-8----
PhenoCode <- read.csv(
    "../UKB/IS-EWAS&MR/data_archived/FinnGen/PhenoCode_5e-8_20241219.csv"
) %>%
    filter(PhenoCode != "") %>%
    dplyr::pull(PhenoCode) %>%
    as.character()
input_dir <- "../UKB/IS-EWAS&MR/data_archived/FinnGen/finngen_R12_5e-8"

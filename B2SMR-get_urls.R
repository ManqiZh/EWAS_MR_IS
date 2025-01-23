# IEU----
get_urls_ieu <- function(ids, out_prefix = "ieu") {
    download_links <- sapply(ids, function(temp_id) {
        file.path(
            "https://gwas.mrcieu.ac.uk/files",
            temp_id, paste0(temp_id, ".vcf.gz")
        )
    })
    fn <- paste0(out_prefix, "_Download_urls.txt")
    write.table(
        download_links,
        file = fn,
        row.names = FALSE,
        quote = FALSE,
        col.names = FALSE
    )
    cli::cli_alert_success(
        "IEU database batch download links output to file: {fn}"
    )
}


# GWAS Catalog----
get_urls_gwas_catalog <- function(data,
                                  download_suffix = "_buildGRCh37.tsv.gz",
                                  out_prefix = "gwas_catalog") {
    download_urls <- sapply(data$summaryStatistics, function(i) {
        temp_list <- strsplit(i, "/")[[1]]
        file.path(i, paste0(temp_list[length(temp_list)], download_suffix))
    })
    fn <- paste0(out_prefix, "_Download_urls.txt")
    write.table(
        download_urls,
        file = fn,
        row.names = FALSE,
        quote = FALSE,
        col.names = FALSE,
        sep = "\t"
    )
    cli::cli_alert_success(
        "GWAS Catalog database batch download links output to file: {fn}"
    )
}


# Finngen----
get_urls_finngen <- function(data, out_prefix = "finngen") {
    download_links <- data$path_https
    fn <- paste0(out_prefix, "_Download_urls.txt")
    write.table(
        download_links,
        file = fn,
        row.names = FALSE,
        quote = FALSE,
        col.names = FALSE
    )
    cli::cli_alert_success(
        "Finngen database batch download links output to file: {fn}"
    )
}


# UKB----
get_urls_ukb <- function(data, out_prefix = "ukb") {
    download_links <- data$AWS.File
    fn <- paste0(out_prefix, "_Download_urls.txt")
    write.table(
        download_links,
        file = fn,
        row.names = FALSE,
        quote = FALSE,
        col.names = FALSE
    )
    cli::cli_alert_success(
        "UKB database batch download links output to file: {fn}"
    )
}

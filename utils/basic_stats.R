get_aroma_stats_by_volatles <- function(gcms_pheno_tbl) {
    tot_sample_ubiq_volatiles <- colSums(gcms_pheno_tbl != 0)
    tot_sample_abund_volatiles <- colSums(gcms_pheno_tbl)
    return(
        data.frame(
            Name      = names(tot_sample_ubiq_volatiles),
            Ubiquity  = tot_sample_ubiq_volatiles,
            Abundance = tot_sample_abund_volatiles
        )
    )
}

get_aroma_stats_by_samples <- function(gcms_pheno_tbl) {
    tot_volatile_ubiq_sample <- rowSums(gcms_pheno_tbl != 0)
    tot_volatile_abund_sample <- rowSums(gcms_pheno_tbl)

    return(
        data.frame(
            Sample      = seq_along(tot_volatile_ubiq_sample),
            Ubiquity  = tot_volatile_ubiq_sample,
            Abundance = tot_volatile_abund_sample
        )
    )
}

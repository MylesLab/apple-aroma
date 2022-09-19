get_aroma_basic_stats_df <- function(gcms_pheno_tbl) {
    #' @description This function is responsible for generating the basic stats
    #'  from GCMS dataframe (without appleid column)
    #' @param gcms_pheno_tbl - GCMS dataframe (without appleid column)
    #' @return a dataframe containing the basic stats


    # get the basic stats
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

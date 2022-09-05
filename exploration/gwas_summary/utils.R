
#' Updates the given top_hits dataframe based on the hits it finds for the compounds in the files
#'
#' @param compounds_list the list of compounds to search for
#' @param files_df the list of files to search those compounds in
#' @param top_hit_df pass-by-value top_hit_df dataframe to accumulate
#' @param compound_peak_map pass-by-value compound_peak_map hash to accumulate
#'
#' @return a data frame containing the top hits for the given compounds
update_top_hits <- function(compounds_list, files_df, top_hits_df, compound_peak_map){
  for(compound_name in compounds_list) {
    fpath <- files_df[grep(compound_name,files_df)]
    print(fpath)
    fs <- fread(fpath)
    pval <- fs$pval

    row <- fs[which(fs$pval == min(fs$pval)),] %>%
      separate(SNP,into = c("CHR","POS","Tool","Allele"), sep = "_", remove = T,extra = "merge") %>%
      mutate(CHR=str_sub(CHR,2,-1)) %>%
      mutate(CHR=as.integer(CHR), POS=as.integer(POS)) %>%
      mutate(CHR=replace(CHR,CHR == 0, 18)) %>%
      mutate(Compound=compound_name)

    # add the compound to the hash
    key <- paste0(row$CHR,"_",row$POS)
    compound_peak_map[[key]] <- unique(c(compound_peak_map[[key]], compound_name))

    # accumulate the top_hits_df
    top_hits_df <- rbind(top_hits_df,row)

  }
  return(
    list(top_hits_df,compound_peak_map)
  )
}
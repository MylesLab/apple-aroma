# Title     : Utility Functions
# Objective : This script contains all the functions used by this module
# Created by: tayabsoomro
# Created on: 2021-07-22

#' Cleanup Compound Names
#'
#' @param df the dataframe whose column names to cleanup 
#'
#' @return a vector of length n (where n is the number of columns in df) which
#'         contains the cleaned names.
#' @export
#'
cleanup_compound_names <- function(df){
  # cleaning up the column names
  new_col_names = NULL
  for(i in colnames(df)){
    new_name <- 
      gsub("([^\\.])\\.([^\\.])\\.{1}([^\\.])","\\1-\\2-\\3",
           gsub("([^\\.])\\.([^\\.])","\\1-\\2",
                gsub("([^\\.])\\.\\.([^\\.])","\\1, \\2",
                     gsub("([^\\.])\\.\\.\\.([^\\.])","\\1, -\\2",
                          gsub("([^\\.])\\.([^\\.])","\\1-\\2",
                               gsub("([^\\.])\\.$","\\1-",
                                    gsub("([^\\.])\\.\\.$","\\1,",
                                         gsub("^X([0-9])","\\1",i)
                                    )
                               )
                          )
                     )
                )
           )
      )
    new_col_names <- c(new_col_names, new_name)
  }
  
  # manually fixing some of the names
  new_col_names[which(new_col_names == "3-Hexenyl-acetate")] <- '3-Hexenyl acetate'
  new_col_names[which(new_col_names == "3-Methyl-3-butenyl-acetate")] <- '3-Methyl-3-butenyl acetate'
  new_col_names[which(new_col_names == "5-ethyl-2-5H, furanone")] <- '5-ethyl-2(5H)-furanone'
  new_col_names[which(new_col_names == "5-Hexenyl-acetate")] <- '5-Hexenyl acetate'
  new_col_names[which(new_col_names == "à-Farnesene")] <- 'a-Farnesene'
  new_col_names[which("X-E, Geranylacetone" == new_col_names)] <- '(E)-Geranylacetone'
  
  # modify the row names with updated names
  return(new_col_names)
  
}
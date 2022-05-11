library(readxl)

# read the compound classificaiton pivot table if not alreay available
if (!exists("classification_pivot")) {
  classification_pivot <- read_excel('data/processed/Supplementary_Data.xlsx', sheet = 'Compound Class')
}


get_compound_class <- function(compound_name) {
    #' Return the compound class a specific compound belongs to
    #'
    #' @param compound_name - the name of the compound
    #'
    #' @return a string value representing the class of compound_name

  classification_pivot[classification_pivot$`Compound name` == compound_name,]$Classification
}


get_all_compounds_with_class <- function(class_name) {
    #' Returns a list of compound names that have the specified compound class
    #'
    #' @param class_name - the name of the class to filter for
    #'
    #' @return a vector containing the names of compounds that belong to specified class

  get("Compound name",classification_pivot[classification_pivot$Classification == class_name,'Compound name'])
}
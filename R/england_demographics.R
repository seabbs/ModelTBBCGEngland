#'  Population data stratified by age group and UK birth status; Extracted, Cleaned, and Aggregated from the LFS
#'
#' A dataset containing population data stratified by birth status, by year from 2000 to 2016, estimmated using the Labour Force Survey.
#'See [tbinenglanddataclean](https://www.samabbott.co.uk/tbinenglanddataclean/).
#' for information regarding the raw data sources used, extraction, and cleaning.
#' @format A data frame with 386 rows and 6 variables.
#' \describe{
#'   \item{year}{Year}
#'   \item{age_group}{Age group to be used in modelling}
#'   \item{population}{Population as estimated from the LFS}
#' }
"england_demographics"

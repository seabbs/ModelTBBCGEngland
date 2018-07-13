#' Yearly Age Specific Mortality Rates
#'
#' A dataset containing estimated age specific mortality rates from 2000 until 2014. These
#' are estimated by the ONS and are 3 year rolling estimates. See [tbinenglanddataclean](www.samabbott.co.uk/tbinenglanddataclean)
#' for details.
#' @format A data frame with - rows and 4 variables.
#' \describe{
#'   \item{year}{Calendar year}
#'   \item{age_group}{5 year age groups (0-4, 5-9, ...., 85-89, 90+)}
#'   \item{mortality_rate}{Estimated age specific mortality rate per year}
#'   \item{exp_life_span}{Estimated expected lifespan if an individual continued to experiance
#'   the same level of age specific mortality}
#' }
"mortality_rates"

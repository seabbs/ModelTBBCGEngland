% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setup_model_obs.R
\name{setup_model_obs}
\alias{setup_model_obs}
\title{Setup Observed Model Data}
\usage{
setup_model_obs(years_of_data = NULL, years_of_age = NULL,
  age_groups = NULL, spacing_of_historic_tb = 1,
  con_age_groups = NULL, aggregated = FALSE, historic = FALSE)
}
\arguments{
\item{years_of_data}{Numeric, the years modern case data to filter for. If not given all are returned.}

\item{years_of_age}{Numeric, the years of age distributed cases to fit to.
If not specified then no years are used.}

\item{age_groups}{Numeric, the numeric age groups to include in the observed data, defaults to \code{NULL} in
which case no age groups are included.}

\item{spacing_of_historic_tb}{Numeric, defaults to 1. Mod to use to identify years of data to use for.}

\item{aggregated}{Logical, defaults to \code{FALSE}. Should aggregated data be used.}

\item{historic}{Logical, defaults to \code{FALSE}. Should historic data be used.}

\item{con_age_age_groups}{Character string, age groups to include as observations. By default no age groups are included. Options
include \code{"children"}, \code{"adults"}, and \code{"older adults"}.}
}
\value{
A named list of observed data required by the model.
}
\description{
Performs data preprocessing required to prepare the observed data
for the model.
}
\examples{

## Code
setup_model_obs

## Output
setup_model_obs(years_of_age = 2000:2015, con_age_groups = c("children", "adults", "older adults"))

}

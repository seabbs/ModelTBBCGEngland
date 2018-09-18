#' Setup Observed Model Data
#'
#' @description Performs data preprocessing required to prepare the observed data 
#' for the model.
#' @param years_of_age Numeric, the years of age distributed cases to fit to. 
#' If not specified then no years are used.
#' @return A named list of observed data required by the model.
#' @export
#'
#' @inheritParams setup_model_input
#' @importFrom dplyr filter group_by mutate summarise select arrange count rename ungroup
#' @examples
#' 
#' ## Code
#' setup_model_obs
#' 
#' ## Output
#' setup_model_obs()
#' 
setup_model_obs <- function(years_of_age = NULL) {
  

  ## Extract historic Pulmonary TB cases
  historic_p_tb <- ModelTBBCGEngland::historic_cases %>%
    filter(year < 2000, year >= 1982) %>% 
    select(time = year, value = pulmonary) %>% 
    mutate(time = time - 1931)
  
  ## Extract age stratified UK born cases
  age_cases <- ModelTBBCGEngland::incidence %>% 
    filter(ukborn == "UK Born") %>% 
    select(-ukborn) %>% 
    group_by(year, age_group) %>% 
    summarise(value = sum(incidence, na.rm = T)) %>% 
    ungroup %>% 
    mutate(time = year - 1931) %>% 
    arrange(time, age_group) %>% 
    mutate(age = as.numeric(age_group) - 1) %>% 
    select(time, age, value) %>% 
    arrange(time, age)
  
  ## Extract UK born cases
  yearly_cases <- age_cases %>% 
    group_by(time) %>% 
    summarise(value = sum(value, na.rm = TRUE))
  
  ## Filter age cases
  if (is.null(years_of_age))  {
    age_cases <- age_cases %>% 
      filter(time %in% (years_of_age - 1931))
  }
  
  obs <- list(
    "YearlyHistPInc" = historic_p_tb,
    "YearlyAgeInc" = age_cases,
    "YearlyInc" = yearly_cases
  )
  
  return(obs)
}

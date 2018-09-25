#' Setup Observed Model Data
#'
#' @description Performs data preprocessing required to prepare the observed data 
#' for the model.
#' @param years_of_data Numeric, the years modern case data to filter for. If not given all are returned.
#' @param years_of_age Numeric, the years of age distributed cases to fit to. 
#' If not specified then no years are used.
#' @param age_groups Numeric, the numeric age groups to include in the observed data, defaults to \code{NULL} in 
#' which case no age groups are included.
#' @param con_age_age_groups Character string, age groups to include as observations. By default no age groups are included. Options 
#' include \code{"children"}, \code{"adults"}, and \code{"older adults"}.
#' @param spacing_of_historic_tb Numeric, defaults to 1. Mod to use to identify years of data to use for.
#' @return A named list of observed data required by the model.
#' @export
#'
#' @inheritParams setup_model_input
#' @importFrom dplyr filter group_by mutate summarise select arrange count rename ungroup add_row
#' @examples
#' 
#' ## Code
#' setup_model_obs
#' 
#' ## Output
#' setup_model_obs(years_of_age = 2000:2015, con_age_groups = c("children", "adults", "older adults"))
#' 
setup_model_obs <- function(years_of_data = NULL,
                            years_of_age = NULL, age_groups = NULL, 
                            spacing_of_historic_tb = 1, 
                            con_age_groups = NULL) {
  

  ## Extract historic Pulmonary TB cases
  historic_p_tb <- ModelTBBCGEngland::historic_cases %>%
    filter(year < 2000, year >= 1982) %>% 
    select(time = year, value = pulmonary) %>% 
    mutate(time = time - 1931) %>% 
    filter((time - min(time)) %% spacing_of_historic_tb == 0)
  
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
  
  if (!is.null(years_of_data)) {
    yearly_cases <- yearly_cases %>% 
      filter(time %in% (years_of_data - 1931))
  }
  
  obs <- list(
    "YearlyHistPInc" = historic_p_tb,
    "YearlyInc" = yearly_cases
  )
  
  ## Filter age cases
  if (!is.null(years_of_age))  {
    
    age_cases <- age_cases %>% 
      filter(time %in% (years_of_age - 1931))
    
    if (!is.null(con_age_groups)) {
      if ("children" %in% con_age_groups) {
        YearlyChildInc <- age_cases %>% 
          filter(age %in% c(0:2)) %>% 
          group_by(time) %>% 
          summarise(value = sum(value))
        
        obs[["YearlyChildInc"]] <- YearlyChildInc
      }
      
      if ("adults" %in% con_age_groups) {
        YearlyAdultInc <- age_cases %>% 
          filter(age %in% c(3:10)) %>% 
          group_by(time) %>% 
          summarise(value = sum(value))
        
        obs[["YearlyAdultInc"]] <- YearlyAdultInc
      }
      
      if ("older adults" %in% con_age_groups) {
        YearlyOlderAdultInc <- age_cases %>% 
          filter(age %in% c(11)) %>% 
          group_by(time) %>% 
          summarise(value = sum(value))
        
        obs[["YearlyOlderAdultInc"]] <- YearlyOlderAdultInc
      }
      
    }
    if (!is.null(age_groups)) {
      age_cases <- age_cases %>% 
        mutate(value = ifelse(!(age %in% age_groups), NA, value))
      
      
      obs[["YearlyAgeInc"]] <- age_cases
    }
    

  }
  
  return(obs)
}

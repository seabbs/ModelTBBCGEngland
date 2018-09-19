#' Setup Input Model Data
#' 
#' @description Performs data preprocessing required to prepare the input data 
#' for the model.
#' @param run_time Numeric, the run time of the model
#' @param time_scale_numeric Numeric the time scale to use (with 1 being a year, 12 a month etc.). Defaults to 1.
#' @return A named list of data inputs required by the model.
#' @export
#'
#' @importFrom dplyr filter group_by mutate summarise select arrange count rename ungroup bind_cols bind_rows n
#' @importFrom tidyr unnest nest
#' @importFrom tibble tibble
#' @importFrom purrr map map2
#' @examples
#' 
#' ## Code
#' setup_model_input
#' 
#' ## Output
#' setup_model_input(run_time = 80, time_scale_numeric = 1)
setup_model_input <- function(run_time = NULL, time_scale_numeric = 1) {
  
  ## Set up initial population distribution
  pop_dist <- england_demographics %>% 
    filter(CoB == "UK born") %>% 
    group_by(year) %>% 
    mutate(age = 0:(n() - 1)) %>% 
    group_by(age) %>% 
    summarise(value = mean(proportion_age_by_year)) %>% 
    select(age, value)
  
  ## Set up births scaling for time horizon
  t_births <- births %>% 
    filter(year >= 1931) %>% 
    mutate(time = year - 1931) %>% 
    mutate(time_n = map(time, ~ tibble(time_n = time_scale_numeric * . + 0:(time_scale_numeric - 1)))) %>% 
    unnest() %>% 
    mutate(value = births / time_scale_numeric) %>% 
    select(time = time_n, value) %>% 
    filter(time <= run_time * time_scale_numeric)
  
  ## Set up expected lifespan
  exp_life_span <- mortality_rates %>% 
    mutate(time = year - 1931,
           value = exp_life_span * time_scale_numeric) %>% 
    group_by(time) %>% 
    mutate(age = 0:(n() - 1)) %>% 
    ungroup %>% 
    mutate(time_n = map(time, ~ tibble(time_n = time_scale_numeric * . + 0:(time_scale_numeric - 1)))) %>% 
    unnest() %>% 
    select(time = time_n, age, value) %>% 
    filter(time <= run_time * time_scale_numeric)
  
  ## Set up Polymod contacts 
  polymod <- contact %>% 
    arrange(age_x, age_y) %>% 
    group_by(age_x) %>% 
    mutate(age2 = 0:(n() - 1)) %>% 
    group_by(age_y) %>% 
    mutate(age = 0:(n() - 1)) %>% 
    ungroup %>% 
    mutate_at(.vars = vars(mean, sd),
              ~ . / time_scale_numeric) %>% 
    select(age, age2, mean, sd)
  
  ## Mean contacts
  polymod_mean <- polymod %>% 
    select(age, age2, value = mean)
  
  ## SD of contacts
  polymod_sd <- polymod %>% 
    select(age, age2, value = sd)
  
  
  ## Extact non-UK born pulmonary cases - estimate previous cases in the model
  nonukborn_p_cases <- incidence %>% 
    filter(ukborn == "Non-UK Born",
           pulmextrapulm == "Pulmonary, with or without EP") %>% 
    select(-ukborn, -pulmextrapulm, -type, -policy_change) %>% 
    mutate(time = year - 1931) %>% 
    arrange(time, age_group) %>% 
    mutate(age = as.numeric(age_group) - 1) %>% 
    select(time, age, incidence) %>% 
    count(time, age, wt = incidence) %>% 
    rename(value = n) %>% 
    mutate(time_n = map(time, ~ tibble(time_n = time_scale_numeric * . + 0:(time_scale_numeric - 1)))) %>% 
    unnest() %>% 
    mutate(value = value / time_scale_numeric) %>% 
    select(time = time_n, age, value)
  
  ## Estimated future non UK born cases using a poisson regression model adjusting for time and age.
    nonukborn_p_cases <-  nonukborn_p_cases  %>% 
    dplyr::filter(time >= (2010 - 1931)) %>% 
    nest() %>% 
    mutate(model = map(data, ~ glm(value ~ time + factor(age),, family = poisson, data = .))) %>% 
    mutate(new_data = map(data,  ~ expand.grid(age = min(.$age):max(.$age), time = c(max(.$time) + 1):(2100 - 1931)))) %>% 
    mutate(pred_cases = map2(model, new_data, ~ tibble(value = predict(.x, .y, type = "response")))) %>% 
    mutate(pred_data = map2(new_data, pred_cases, ~ bind_cols(.x, .y))) %>% 
    mutate(all_data = map(pred_data, ~ bind_rows(nonukborn_p_cases, .))) %>% 
    select(all_data) %>% 
    unnest(all_data)
  
  
  ## Non UK born cases in 2000 - used to estimate historic non UK born cases
  NUKCases2000 <- nonukborn_p_cases %>% 
    filter(time == time_scale_numeric * (2000 - 1931)) %>% 
    select(-time)
  
  
  
  input <- list(
    "pop_dist" = pop_dist,
    "births_input" = t_births,
    "exp_life_span" = exp_life_span,
    "polymod" = polymod_mean,
    "polymod_sd" = polymod_sd,
    "NonUKBornPCases" = nonukborn_p_cases,
    "NUKCases2000" = NUKCases2000
  )  
  
  return(input)
}
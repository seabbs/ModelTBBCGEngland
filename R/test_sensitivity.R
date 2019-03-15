#' Test the sensitivity of a Libbi model using PRCC
#'
#' @param model A LiBBi model object
#' @param obs A character string containing the observation to use to calculate the PRCC.
#' @param target_time Numeric, the time at which to estimate the model sensitivity. If not specified
#' then this will default to the last fitted point that the model has produced output for.
#'
#' @return A data frame containing the names of the parameters in the model, the correlation with the outcome
#' and the p value of this correlation
#' @export
#' @importFrom rbi bi_read
#' @importFrom dplyr bind_rows select arrange
#' @importFrom tidyr spread
#' @importFrom epiR epi.prcc
#' @examples
#' 
#' test_sensitivity
test_sensitivity <- function(model = NULL, obs = NULL, target_time = NULL) {
  
  model_params <- rbi::bi_read(model, type = "param") %>% 
    bind_rows(.id = "Parameter")
  
  model_obs <- rbi::bi_read(model) %>% 
    {.[names(.) %in% obs]} %>%  
    bind_rows(.id = "Observation") 
  if (is.null(target_time)) {
    target_time <- model_obs$time %>% 
      max
  }
  
  obs <- model_obs %>% 
    dplyr::filter(time == target_time)

  
  
  ## Spread parameters into the form the prcc funciton expects
  params <- model_params %>% 
    select(Parameter, value, np) %>% 
    spread(key = "Parameter", value = "value") %>% 
    select(-np)
  
  ## Join obs to parameters
  sample <- params %>% 
    bind_cols(obs %>% 
                select(value) %>% 
                setNames("Observation"))
  
  prcc <- epi.prcc(sample) %>% 
    mutate(Parameter = colnames(params)) %>% 
    select(Parameter, gamma, p.value) %>% 
    arrange(desc(abs(gamma)))
  
  
  return(prcc)
}

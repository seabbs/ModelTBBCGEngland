#' Test Model
#'
#' @description A function to test the models in development through the full fitting pipeline based on synthetic data.
#' @param model A character string containing the model name. Alternatively a model loaded via `rbi::bi_model` may be passed.
#' @param gen_data Logical, defaults to /code{TRUE}. Should data be synthesised using the model and priors.
#' @param run_time Numeric, the number of years to run the model fitting and simulation for, defaults to 74 (i.e from 1931 until 2005).
#' @param time_scale Character, defaults to \code{"year"}. A monthly timescale can also be set with \code{"month" }.
#' @param plot_input_data Logical, defaults to \code{TRUE}. Should input data be plotted
#' @param sample_priors Logical, defaults to \code{FALSE}. Should the model priors be sampled.
#' @param prior_samples Numeric, the number of samples to take from the priors. Defaults to 1000.
#' @param posterior_samples Numeric, the number of samples to take from the posterior estimated using pmcmc (requires \code{fit = TRUE}). Defaults to 1000.
#' @param nparticles Numeric, the initial number of particles to use in the particle filters. Defaults to \code{NULL}, in which case the number of outputs 
#' is used as the initial particle number.
#' @param nthreads Numeric, defaults to 4. The number of parallel jobs to run. The most efficient option is likely to be to match the 
#' number of cores available.
#' @param adapt_particles Logical, defaults to \code{FALSE}. Should the number of particles be adapted.
#' @param adapt_proposal  Logical, defaults to \code{TRUE}. Should the proposal be adjusted to improve the acceptance rate.
#' @param min_acc Numeric, defaults to 0.05. The minimum target acceptance rate.
#' @param max_acc Numeric, defaults to 0.4. The maximum target acceptance rate.
#' @param fit Logical, defaults to \code{TRUE}. Should the model be fitted with 1000 samples extracted.
#' @param save_output Logical, defaults to \code{FALSE}. Should the model results be saved as a test dataset.
#' @param verbose Logical, defaults to \code{TRUE}. Should progress messages and output be printed.
#' @param libbi_verbose Logical, defaults to \code{FALSE}. Should \code{libbi} produce verbose progress and warnings messages.
#' @param browse Logical, defaults to \code{FALSE}. Should the function be run in debug mode using \code{browser}.
#' @param const_pop Logical, defaults to \code{FALSE}. Should a constant population be maintained using births equals deaths.
#' @param no_age Logical, defaults to \code{FALSE}. Should ageing be disabled.
#' @param no_disease Logical, defaults to \code{FALSE}. Should disease be removed from the model
#' @return A LibBi model object based on the inputed test model.
#' @export
#'
#' @importFrom rbi fix bi_model sample bi_read bi_generate_dataset libbi get_block save_libbi
#' @import rbi.helpers 
#' @import ggplot2
#' @importFrom dplyr filter mutate select vars arrange
#' @importFrom stats runif time
#' @importFrom utils str
#' @importFrom graphics plot
#' @importFrom purrr map
#' @importFrom tibble tibble
#' @importFrom tidyr unnest
#' @examples
#' 
#' 
test_model <- function(model= "BaseLineModel", gen_data = TRUE, run_time = 74, time_scale = "year", plot_input_data = TRUE,
                       sample_priors = TRUE, prior_samples = 1000, nparticles = NULL, adapt_particles = FALSE,
                       adapt_proposal = FALSE, min_acc = 0.05, max_acc = 0.4, fit = FALSE, posterior_samples = 1000, 
                       save_output = FALSE, nthreads = 4, verbose = TRUE, libbi_verbose = FALSE, browse = FALSE,
                       const_pop = FALSE, no_age = FALSE, no_disease = FALSE) {
  
  

  # Set the time scale for the model ----------------------------------------


  if (time_scale == "year") {
    time_scale_numeric <- 1
  }else if (time_scale == "month") {
    time_scale_numeric <- 12
  }else{
    stop("Only a yearly (year) or monthly (month) timescale have been implemented.")
  }
  

  # Set the number of particles ---------------------------------------------

  if (is.null(nparticles)) {
    nparticles <- run_time * time_scale_numeric
  }  
  
  # Load the model ----------------------------------------------------------
  
  if (is.character(model)) {
    model_file <- system.file(package="ModelTBBCGEngland", paste0("bi/", model, ".bi"))
    
    tb_model_raw <- bi_model(model_file)
    
  }else{
    tb_model_raw <- model
  }
  
  

# Specify model setup conditions ------------------------------------------

  if (time_scale == "year") {
    tb_model_raw <- fix(tb_model_raw, ScaleTime = 1 / 12)
  }else if (time_scale == "month") {
    tb_model_raw <- fix(tb_model_raw, ScaleTime = 1)
  }else{
    stop("Only a yearly (year) or monthly (month) timescale have been implemented.")
  }
  
  if (const_pop) {
    tb_model_raw <- fix(tb_model_raw, const_pop = 1)
  }
  
  if (no_age) {
    tb_model_raw <- fix(tb_model_raw, no_age = 1)
  }
  
  if (no_disease) {
    tb_model_raw <- fix(tb_model_raw, no_disease = 1)
  }
  
  # Allow for interactive debugging -----------------------------------------
  
  if (browse) {
    browser() 
  }


# Set up model input ------------------------------------------------------
  
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
      rename(value = n)
  
  ## Estimate non-ukborn cases from 1982 until 1999 - do in model
      
    
input <- list(
  "pop_dist" = pop_dist,
  "births_input" = t_births,
  "exp_life_span" = exp_life_span,
  "polymod" = polymod_mean,
  "polymod_sd" = polymod_sd,
  "NonUKBornPCases" = nonukborn_p_cases
)  


# Set up abs data ---------------------------------------------------------

  ## Extract historic Pulmonary TB cases
  historic_p_tb <- historic_cases %>%
    filter(year < 2000) %>% 
    select(time = year, value = pulmonary) %>% 
    mutate(time = time - 1931)
  
  ## Extract age stratified UK born cases
  age_cases <- incidence %>% 
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

  obs <- list(
    "YearlyHistPInc" = historic_p_tb,
    "YearlyAgeInc" = age_cases,
    "YearlyInc" = yearly_cases
  )
  # Generate data from the model --------------------------------------------
  if (gen_data) {
    
    if (verbose) {
      message("Generating data from model")
    }
    
    tb_data <- bi_generate_dataset(tb_model_raw, end_time = run_time * time_scale_numeric, 
                                   input = input, obs = obs, 
                                   noutputs =  run_time * time_scale_numeric,
                                   seed = 1234, verbose = libbi_verbose)
    
    if (verbose) {
      message("Summary of generated model data")
      print(tb_data)
    }
    
    
    obs <- bi_read(tb_data) %>% 
      .[names(obs)]
    
    if (verbose) {
      message("Variables in the generated data")
      print(names(obs))
    }
    
  }
  
  
  # Plot incidence generated by the model -----------------------------------
  if (plot_obs) {
    
    time <- NULL; value <- NULL;
    
    message("Plot historic number of pulmonary cases detected each year:")
    p_hist_cases <- obs$YearlyHistPInc %>% 
      dplyr::filter(time > 0) %>% 
      ggplot(aes(x = time, y = value)) +
      geom_point(size = 1.2) +
      geom_line(size = 1.1, alpha = 0.6) +
      theme_minimal() +
      labs(x = "Time",
           y = "Yearly notified cases")
    
    print(p_hist_cases)
    
    
    message("Plot number of cases detected each year:")
    p_cases <- obs$YearlyCases %>% 
      dplyr::filter(time > 0) %>% 
      ggplot(aes(x = time, y = value)) +
      geom_point(size = 1.2) +
      geom_line(size = 1.1, alpha = 0.6) +
      theme_minimal() +
      labs(x = "Time",
           y = "Yearly notified cases")
    
    print(p_cases)
    
    
    message("Plot number of cases detected each year, stratified by age group:")
    p_age_cases <- obs$YearlyAgeCases %>% 
      dplyr::filter(time > 0) %>% 
      ggplot(aes(x = time, y = value)) +
      geom_point(size = 1.2) +
      geom_line(size = 1.1, alpha = 0.6) +
      theme_minimal() +
      labs(x = "Time",
           y = "Yearly notified cases") +
      facet_wrap(~age)
    
    print(p_age_cases)
    
  }

  
  # Load the model ----------------------------------------------------------
  
  if (verbose) {
    message("Load the model as a compiled libbi object")
  }

  tb_model <- libbi(tb_model_raw, input = input, obs = obs)
  
  
  # Sample from priors ------------------------------------------------------
if (sample_priors) {
  if (verbose) {
    message("Sample priors")
  }
  
  tb_model <- sample(tb_model, target="prior", nsamples = prior_samples, 
                     end_time = run_time* time_scale_numeric, 
                     noutputs = run_time* time_scale_numeric,
                     nparticles = nparticles, nthreads = nthreads, verbose = libbi_verbose)
  
  if (verbose) {
    message("Summary of prior sampling")
    print(tb_model)
    message("Load priors")
  }
  
  priors <- bi_read(tb_model)
  
  if (verbose) {
    message("Structure of priors")
    str(priors)
    
    message("Plot priors")
    print(plot(tb_model))
    
  }
  
  
}  

  # Adapting the number of particles ----------------------------------------
  
  if (adapt_particles) {
    if (verbose) {
      message("Adapting particles")
    }
    
    tb_model <- adapt_particles(tb_model)
  }
  
  
  # Adapting the proposal ---------------------------------------------------
  
  if (adapt_proposal) {
    if (verbose) {
      message("Adapting particles with a min acceptance of ", min_acc, " and a maximum acceptance of ", max_acc)
    }
    
    tb_model <- adapt_proposal(tb_model, min = min_acc, max = max_acc)
    
    if (verbose) {
      message("Adapted proposal distribution")
      prop_block <- get_block(tb_model$model, "proposal_parameter")
      
      print(prop_block)
    }
    
  }
  
  
  # Test model fitting ------------------------------------------------------
  
  
  if (fit) {
    if (verbose) {
      message("Fitting using PMCMC")
    }
    
    
    tb_model <- tb_model %>% 
      sample(target="posterior", sample_obs = TRUE, nsamples = posterior_samples,
             noutputs = run_time* time_scale_numeric,
             verbose = libbi_verbose, nthreads = nthreads)
    
    
    
    # Look at chain -----------------------------------------------------------
    
    p <- plot(tb_model, plot = FALSE)
    
    if (verbose) {
      
      message("Plotting Trajectories")
      print(p$trajectories)
      
      message("Plotting traces")
      print(p$traces)
      
      message("Plotting densities")
      print(p$densities)
      
      message("Overview of all data")
      print(p$data)
    }
  }
  
  
  
  if (save_output) {
    if (verbose) {
      message("Saving model output")
    }
    
    test_data_loc <- paste0("vignettes/results/model-", Sys.Date(), "-", round(runif(1, 0, 10000)), ".rds")
    save_libbi(tb_model, test_data_loc)
    
    if (verbose) {
      message("Test data saved to: ")
      message(test_data_loc)
    }
  }
  return(tb_model)
}
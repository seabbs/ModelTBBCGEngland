## This script is used to fit all scenarios considered including the baseline scenario
## Options exist to specify which scenarios to consider and whether to fit them sequentially or in parallel.

# Analysis settings -------------------------------------------------------

cores <- future::availableCores()[[1]] ## Cores to use for analysis, defaults to all detected.
parallel_scenarios <- 1 ##Number of scenarios to fit in parallel. If set to be higher than 1 then 
                        ##each scenario uses cores / parallel_scenarios (rounded down).
scenario <- "baseline"   ##Named scenario to evaluate.
dir_path <- "./vignettes/results" ##Path to results, these folders must exist and be writable.
calib_run <-  FALSE
sample_priors <- FALSE
adapt_part <- FALSE
adapt_prop <- FALSE
fit <- FALSE

GetoptLong::GetoptLong(
  "cores=f", "Number of cores to use for evaluation, defaults to all detected cores.",
  "parallel_scenarios=f", "Number of scenarios to run in parallel, defaults to 1.", 
  "scenario=s@", "Named scenarios to evaluate pass multiple scenarios using this arguement each time. Defaults to all scenarios",
  "dir_path=s", "Directory to save the evaluated scenarios into. Defaults to ./vignettes/results",
  "calib_run", "Should the model be run for calibration. This uses a reduced set of data points. The results may be used to inform the main model runs.",
  "sample_priors", "Should priors be sampled",
  "adapt_part", "Should the number of particles be adapted",
  "adapt_prop", "Should the proposal distribution be adapted",
  "fit", "Should the scenarios be fitted"
)


# Set calibration options -------------------------------------------------

## Give calibration run settings
if (calib_run) {
  years_of_data <- 2000
  years_of_age <- NULL
  nparticles <- cores
  run_time <- 69
  adapt_part_samples <- 100
  adapt_prop_samples <- 100
  adapt_part_it <- 5
  adapt_prop_it <- 3
  adapt_scale <- 2
}else{
  years_of_data <- NULL
  years_of_age <- c(2000, 2004)
  nparticles <- NULL
  run_time <- 73
  adapt_part_samples <- 100
  adapt_prop_samples <- 100
  adapt_part_it <- 3
  adapt_prop_it <- 2
  adapt_scale <- 1.5
}

# Set up analysis ---------------------------------------------------------

scenario_path <- paste0(dir_path, "/evaluated-scenarios", "-", stringr::str_replace(Sys.time(), " ", "_"))

if (calib_run) {
  scenario_path <- paste0(scenario_path, "-calibration-run")
}

if (!dir.exists(scenario_path)) {
  dir.create(scenario_path)
}

message("Scenarios output to: ", scenario_path)

## Make the log file
logs <- file.path(scenario_path, "log.txt")
con <- file(logs)

## Send Output to log
sink(con, append = TRUE)
sink(con, type = "message", append = TRUE)


# Load packages -----------------------------------------------------------

library(ModelTBBCGEngland)
library(rbi.helpers)
library(tidyverse)
library(furrr)


# Set up parallel processing ----------------------------------------------

## Set up the number of cores to use for each process
nthreads <- floor(cores / parallel_scenarios)

## Set up processing plan
if (parallel_scenarios == 1) {
  plan(sequential)
}else{
  plan(multiprocess, workers = parallel_scenarios)
}
# Set up model fitting for all scenarios ----------------------------------
## Arguements required for all scenarios: dir_name, scenario specific settings
fit_model_with_baseline_settings <- partial(fit_model,
                                            ## Run time arguements
                                            model = "BaseLineModel", gen_data = FALSE, run_time = run_time,
                                            time_scale = "year", plot_obs = TRUE, nthreads = nthreads,
                                            ##Prior settings
                                            sample_priors = sample_priors, prior_samples = 1000,
                                            ##Particle settings
                                            adapt_particles = adapt_part, nparticles = nparticles, adapt_part_samples = adapt_part_samples ,
                                            adapt_part_it = adapt_part_it, 
                                            ##Proposal settings
                                            adapt_proposal = adapt_prop, adapt_prop_samples = adapt_prop_samples, adapt_prop_it = adapt_prop_it, 
                                            adapt = "shape", adapt_scale = adapt_scale, min_acc = 0.2, max_acc = 0.4,
                                            ##Posterior sampling settings
                                            fit = fit, posterior_samples = 5000, sample_ess_at = 0.2,
                                            rejuv_moves = NULL,
                                            ##Prediction settings
                                            pred_states = FALSE,
                                            ## Model settings
                                            scale_rate_treat = TRUE,
                                            years_of_data = years_of_data, 
                                            years_of_age = years_of_age, 
                                            age_groups = NULL, con_age_groups = c("children", "older adults"), 
                                            spacing_of_historic_tb = 10, noise = TRUE, 
                                            ##Results handling settings)
                                            verbose = TRUE, libbi_verbose = TRUE, 
                                            fitting_verbose = TRUE, save_output = TRUE, 
                                            dir_path = scenario_path, reports = TRUE)




# Outline Scenarios -------------------------------------------------------
scenarios <- list()

## Baseline scenario: Linear scaling for non-UK born cases, homogeneous non-UK born mixing and constant transmission probability across all age groups
scenarios$baseline <- list(
  dir_name = "baseline"
)


# Non-UK born mixing ------------------------------------------------------


## Hetergeneous non-UK born mixing
scenarios$het_non_uk <- list(
  dir_name = "het_non_uk",
  non_uk_mixing = "het"
)


# Transmission probability degrees of freedom -----------------------------


##Variable transmission probability between children and adults
scenarios$trans_prob_var_children <- list(
  dir_name = "trans_prob_var_children",
  trans_prob_freedom = "child_free"
)

##Variable transmission probability between children, older adults and adults
scenarios$trans_prob_var_children_older_adults <- list(
  dir_name = "trans_prob_var_children_older_adults",
  trans_prob_freedom = "child_older_adult_free"
)


# Non-UK born scaling -----------------------------------------------------

##Log / log(max) scaling of non-UK born cases
scenarios$log_non_uk <- list(
  dir_name = "log_non_uk",
  non_uk_scaling = "log"
)

## Constant non UK born cases
scenarios$constant_non_uk <- list(
  dir_name = "constant_non_uk",
  non_uk_scaling = "constant"
)

##  Filter for selected scenarios.
if (!is.null(scenario)) {
  scenarios <- scenarios[scenario]
}

# Set up scenario evaluation ----------------------------------------------
##Requires a list of optional settings to pass to the fit_model function
##All other options given above
evaluate_scenario <- function(scenario) {
  message("Evaluating ", scenario$name, "at ", Sys.time())
  
  model <- do.call(fit_model_with_baseline_settings, scenario)
  
  ## Evaluate model fit via DIC
  if (fit) {
    dic <- DIC(model)
    
    ## Report model DIC
    message(scenario$name, "DIC: ", dic)
  }else{
    message("Model not fitted and therefore DIC not evaluated")
    dic <- NULL
  }

  
  ## Return model DIC
  return(dic)
  
}


# Fit scenarios -----------------------------------------------------------

fitted_scenarios <- future_map_dfr(scenarios, evaluate_scenario, .id = "scenario")   

message("Scenario evaluation complete")

saveRDS(fitted_scenarios, file.path(scenario_path, "scenario_dics.rds"))

# Wind up script ----------------------------------------------------------

sink(file = NULL) 
sink(file = NULL, type = "message")

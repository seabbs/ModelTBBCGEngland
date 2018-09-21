## This script is used to fit all scenarios considered including the baseline scenario
## Options exist to specify which scenarios to consider and whether to fit them sequentially or in parallel.


# Load packages required --------------------------------------------------
library(ModelTBBCGEngland)
library(rbi.helpers)
library(tidyverse)
library(furrr)


# Analysis settings -------------------------------------------------------

cores <- future::availableCores()[[1]] ## Cores to use for analysis, defaults to all detected.
parallel_scenarios <- 1 ##Number of scenarios to fit in parallel. If set to be higher than 1 then 
                        ##each scenario uses cores / parallel_scenarios (rounded down).
scenario <- NULL   ##Named scenario to evaluate.
dir_path <- "./vignettes/results" ##Path to results, these folders must exist and be writable.


GetoptLong::GetoptLong(
  "cores=f", "Number of cores to use for evaluation, defaults to all detected cores.",
  "parallel_scenarios=f", "Number of scenarios to run in parallel, defaults to 1.", 
  "scenario=s@", "Named scenarios to evaluate pass multiple scenarios using this arguement each time. Defaults to all scenarios",
  "dir_path=s", "Directory to save the evaluated scenarios into. Defaults to ./vignettes/results"
)
# Set up analysis ---------------------------------------------------------

## Set up processing plan
if (parallel_scenarios == 1) {
  plan(sequential)
}else{
  plan(multiprocess(workers = parallel_scenarios))
}
## Set up the number of cores to use for each process
nthreads <- floor(cores / parallel_scenarios)

scenario_path <- paste0("vignettes/results/evaluated-scenarios", "-", str_replace(Sys.time(), " ", "_"))

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

# Set up model fitting for all scenarios ----------------------------------
## Arguements required for all scenarios: dir_name, scenario specific settings

fit_model_with_baseline_settings <- partial(fit_model,
                                            ## Run time arguements
                                            model = "BaseLineModel", gen_data = FALSE, run_time = 73,
                                            time_scale = "year", plot_obs = TRUE, nthreads = nthreads,
                                            ##Prior settings
                                            sample_priors = TRUE, prior_samples = 10,
                                            ##Particle settings
                                            adapt_particles = FALSE, nparticles = NULL, adapt_part_samples = 250,
                                            adapt_part_it = 3, 
                                            ##Proposal settings
                                            adapt_proposal = FALSE, adapt_prop_samples = 250, adapt_prop_it = 4, 
                                            adapt = "size", adapt_scale = 1.2, min_acc = 0.05, max_acc = 0.3,
                                            ##Posterior sampling settings
                                            fit = TRUE, posterior_samples = 10, sample_ess_at = 0.5,
                                            rejuv_moves = 5,
                                            ##Prediction settings
                                            pred_states = TRUE,
                                            ## Model settings
                                            scale_rate_treat = TRUE, years_of_age = c(2000, 2004),
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

## Hetergeneous non-UK born mixing
scenarios$het_non_uk <- list(
  dir_name = "het_non_uk",
  non_uk_mixing = "het"
)

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
  
  model <- do.call(fit_model_with_baseline_settings, scenario)
  
  ## Evaluate model fit via DIC
  dic <- DIC(model)
  
  ## Report model DIC
  message(scenario$name, "DIC: ", dic)
  
  ## Return model DIC
  return(dic)
  
}


# Fit scenarios -----------------------------------------------------------

fitted_scenarios <- future_map_dfr(scenarios, evaluate_scenario, .id = "scenario", .progress = TRUE)   

message("Scenario evaluation complete")

saveRDS(fitted_scenarios, file.path(scenario_path, "scenario_dics.rds"))

# Wind up script ----------------------------------------------------------

sink(file = NULL) 
sink(file = NULL, type = "message")

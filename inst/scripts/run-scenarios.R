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

scenario_path <- paste0("evaluated-scenarios", "-", str_replace(Sys.time(), " ", "_"))

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
                                            sample_priors = TRUE, prior_samples = 1000,
                                            ##Particle settings
                                            nparticles = NULL, adapt_particles = FALSE, adapt_part_samples = 100,
                                            adapt_part_it = 1, 
                                            ##Proposal settings
                                            adapt_proposal = TRUE, adapt_prop_samples = 250, adapt_prop_it = 4, 
                                            adapt = "both", adapt_scale = 2, min_acc = 0.2, max_acc = 0.4,
                                            ##Posterior sampling settings
                                            fit = TRUE, posterior_samples = 5000, sample_ess_at = 0.8,
                                            rejuv_moves = NULL,
                                            ##Prediction settings
                                            pred_states = TRUE,
                                            ## Model settings
                                            scale_rate_treat = TRUE, years_of_age = c(2000, 2004),
                                            noise = TRUE, 
                                            ##Results handling settings)
                                            verbose = TRUE, libbi_verbose = TRUE, 
                                            fitting_verbose = TRUE, save_output = TRUE, 
                                            dir_path = scenario_path, reports = TRUE)




# Outline Scenarios -------------------------------------------------------
scenarios <- list()

scenarios$baseline <- list(
  dir_name = "baseline"
)

##  Filter for selected scenarios.
if (!is.null(scenario)) {
  scenarios <- scenarios[[scenario]]
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

fitted_scenarios <- future_map_dfr(scenarios, evaluate_scenario, .id = "scenario")   

message("Scenario evaluation complete")

saveRDS(fitted_scenarios, file.path(scenario_path, "scenario_dics.rds"))

# Wind up script ----------------------------------------------------------

sink(file = NULL) 
sink(file = NULL, type = "message")

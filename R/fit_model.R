#' Fit Model
#'
#' @description A function to fit the models in development through the full fitting pipeline based on synthetic data or observed data.
#' @param model A character string containing the model name. Alternatively a model loaded via `rbi::bi_model` may be passed.
#' @param previous_model_path A character string giving the path to a previous LiBBi model (saved as an .rds). This will replace the default model and 
#' override initial settings.
#' @param gen_data Logical, defaults to /code{TRUE}. Should data be synthesised using the model and priors.
#' @param run_time Numeric, the number of years to run the model fitting and simulation for, defaults to 74 (i.e from 1931 until 2005).
#' @param time_scale Character, defaults to \code{"year"}. A monthly timescale can also be set with \code{"month" }.
#' @param plot_obs Logical, defaults to \code{TRUE}. Should input data be plotted
#' @param sample_priors Logical, defaults to \code{FALSE}. Should the model priors be sampled.
#' @param prior_samples Numeric, the number of samples to take from the priors. Defaults to 1000.
#' @param posterior_samples Numeric, the number of samples to take from the posterior estimated using pmcmc (requires \code{fit = TRUE}). Defaults to 1000.
#' @param nparticles Numeric, the initial number of particles to use in the particle filters. Defaults to \code{NULL}, in which case the number of data points 
#' is used as the initial particle number rounded to the nearest power of 2.
#' @param nthreads Numeric, defaults to 4. The number of parallel jobs to run. The most efficient option is likely to be to match the 
#' number of cores available.
#' @param optim Logical defaults to \code{TRUE}. Should the determinsitic model be optimised prior to other fitting steps.
#' @param adapt_particles Logical, defaults to \code{FALSE}. Should the number of particles be adapted.
#' @param adapt_part_samples Numeric, the number of samples to take from the priors when adapting the number of particles. Defaults to 1000.
#' @param adapt_part_it Numeric, the number of iterations to use for adapting the number of particles. Defaults to 10.
#' @param min_particles Numeric, the starting number of particles to use for adaption. If \code{NULL} uses \code{nparticles}.
#' @param max_particles Numeric, the maximum number of particles to use for adaption. If \code{NULL} then defaults to \code{4 * min_particles}
#' @param proposal_param_block A character string containing a LiBBi proposal parameter block. Defaults to \code{NULL} in 
#' which case the LiBBi defaults will be used.
#' @param proposal_initial_block A character string containing a LiBBi proposal initial block. Defaults to \code{NULL} in 
#' which case the LiBBi defaults will be used.
#' @param adapt_proposal  Logical, defaults to \code{TRUE}. Should the proposal be adjusted to improve the acceptance rate.
#' @param adapt_prop_samples Numeric, the number of samples to take when adapting the proposal distribution. Defaults to 1000. Initially 5 times this number will be sampled
#' in order to move the chain out of high unlikely parameter regions.
#' @param adapt_prop_it Numeric, defaults to 10. The number of iterations to use for adapting the proposal.
#' @param adapt Character string, defaults to "both". The type of adaption to use for the proposal see \code{rbi.helpers::adapt_proposal} for details.
#' @param adapt_scale Numeric, defaults to 2. The scale to use to adapt the proposal.
#' @param min_acc Numeric, defaults to 0.05. The minimum target acceptance rate.
#' @param max_acc Numeric, defaults to 0.4. The maximum target acceptance rate.
#' @param fit Logical, defaults to \code{TRUE}. Should the model be fitted with 1000 samples extracted.
#' @param thin Numeric, the thinning interval to use between posterior samples. May reduce the correlation between samples and reduces memory.
#' @param burn_prop Numeric (between 0 and 1). The proportion of the chain to discard as burn-in.
#' @param sample_ess_at Numeric defaults to 0.8. The thresold of the effective sample size (ess) at which to rejuvernate the particles.
#' @param rejuv_moves Numeric, defaults to \code{NULL}. The number of MCMC samples to take for each rejuvernation step. If \code{NULL} is set so that
#' the acceptance rate of each rejuvernation step is at least 1 minus the effective sample size thresold. This is based on the estimated acceptance rate after adaption. If proposal adaption is not used then it assumed that
#' the accpetance rate is 5\% only usable for testing purposes.
#' @param save_output Logical, defaults to \code{FALSE}. Should the model results be saved as a test dataset.
#' @param verbose Logical, defaults to \code{TRUE}. Should progress messages and output be printed.
#' @param libbi_verbose Logical, defaults to \code{FALSE}. Should \code{libbi} produce verbose progress and warnings messages.
#' @param fitting_verbose Logical, defaults to \code{FALSE}. Should \code{libbi} produce verbose progress and warnings messages whilst fitting.
#' @param browse Logical, defaults to \code{FALSE}. Should the function be run in debug mode using \code{browser}.
#' @param const_pop Logical, defaults to \code{FALSE}. Should a constant population be maintained using births equals deaths.
#' @param no_age Logical, defaults to \code{FALSE}. Should ageing be disabled.
#' @param no_disease Logical, defaults to \code{FALSE}. Should disease be removed from the model
#' @param scale_rate_treat Logical, defaults to \code{TRUE}. Scales the rate of treatment over time from introduction in 1952 to assumed modern standard in 1990.
#' @param years_of_age Numeric, the years of age distributed cases to fit to. Defaults to all years available.
#' @param noise Logical, should process noise be included. Defaults to \code{TRUE}. If \code{FALSE} then noise will still be included 
#' from the measurement model.
#' @param initial_uncertainty Logical, should initial state and parameter uncertainty be included. Defaults to \code{TRUE}. 
#' @param pred_states Logical defaults to \code{TRUE}. Should states be predicted over all time (from model initialisation to 35 years ahead of final run time). 
#' If set to \code{FALSE} states will only be estimated for times with observed data points.
#' @param scale_transmission A character string defaults to \code{"none"}. The default ensures that the transmission probability is constant across all age groups. Other options include;
#' \code{"young_adult"}. These add a modifiying parameter for young adults (15-24).
#' @param scale_noukborn_mixing  A character string defaults to \code{"none"}. The default ensures that non-UK born mixing is constant across all age groups. Other options include;
#' \code{"young_adult"}. These add a modifiying parameter for young adults (15-24).
#' @param seed Numeric, the seed to use for random number generation.
#' @param reports Logical, defaults to \code{TRUE}. Should model reports be generated. Only enabled when \code{save_output = TRUE}.
#' @param time_for_resampling Numeric, defaults to 0 (i.e off). Overall real time (minutes) to allocate to move steps  for the SMC sampler. If set to be non-zero then this will
#' overvide rejuvernaiton and effective sample size setting.
#' @param aggregated_observed Logical, defaults to \code{FALSE}. Should aggregated observational data be used. 
#' @return A LibBi model object based on the inputed test model.
#' @export
#' @inheritParams setup_model_obs
#' @importFrom rbi fix bi_model sample bi_read bi_generate_dataset libbi get_block optimise sample_obs
#' @import rbi.helpers 
#' @import ggplot2
#' @importFrom dplyr filter mutate select vars arrange count rename summarise
#' @importFrom tidyr drop_na
#' @importFrom stats runif time
#' @importFrom utils str
#' @importFrom graphics plot
#' @importFrom purrr map map_dbl map_dfr
#' @importFrom tibble tibble
#' @importFrom stringr str_replace
#' @importFrom rmarkdown render
#' @examples
#' 
#' ## Function code:
#' fit_model
fit_model <- function(model = "BaseLineModel", previous_model_path = NULL, gen_data = TRUE, 
                      run_time = 73, time_scale = "year", plot_obs = TRUE,
                      sample_priors = TRUE, prior_samples = 1000, optim = TRUE, nparticles = NULL,
                      adapt_particles = TRUE, adapt_part_samples = 1000, adapt_part_it = 10, 
                      min_particles = NULL, max_particles = NULL,
                      proposal_param_block = NULL, proposal_initial_block = NULL, 
                      adapt_proposal = TRUE, adapt_prop_samples = 100, adapt_prop_it = 3, adapt = "both",
                      adapt_scale = 2, min_acc = 0.04, max_acc = 0.4,
                      fit = FALSE, posterior_samples = 10000, thin = 1, burn_prop = 0, time_for_resampling = 0, 
                      sample_ess_at = 0.8,
                      rejuv_moves = NULL, nthreads = parallel::detectCores(), verbose = TRUE, libbi_verbose = FALSE, 
                      fitting_verbose = TRUE, pred_states = TRUE, browse = FALSE,
                      const_pop = FALSE, no_age = FALSE, no_disease = FALSE, scale_rate_treat = TRUE, years_of_data = c(2000:2004),
                      years_of_age = c(2000, 2004), age_groups = NULL, con_age_groups = NULL, spacing_of_historic_tb = 10,
                      aggregated_observed = FALSE, initial_uncertainty = TRUE, noise = TRUE,
                      scale_transmission = "none", scale_nonukborn_mixing = "none",
                      save_output = FALSE, dir_path = NULL, dir_name = NULL, reports = TRUE,
                      seed = 1234) {


# Util functions ----------------------------------------------------------

# Check parameters --------------------------------------------------------

if(burn_prop > 1 || burn_prop < 0) {
  stop("The burn in proportion must be less or equal to 1 and greater than or equal to 0.")
}  

  # Set up model directory --------------------------------------------------
if (save_output) {
  if (is.null(dir_path)) {
    dir_path <- "."
  }
  
  if (is.null(dir_name)) {
    dir_name <- "model-run"
  }
  
  ## Add time to filename
  dir_name <- paste0(dir_name, "-", str_replace(Sys.time(), " ", "_"))
  
  ## Construct path
  model_dir <- file.path(dir_path, dir_name)
  
  message("Model output: ", model_dir)
  
  ## Make the folder
  if (!dir.exists(model_dir)) {
    dir.create(model_dir)
    
    ## Make a sub folder for data
    data_dir <- file.path(model_dir, "data")
    
    dir.create(data_dir)
    
    ## Make a sub folder for plots
    plots_dir <- file.path(model_dir, "plots")
    
    dir.create(plots_dir)
    
    ## Make a sub folder for libbi output
    libbi_dir <- file.path(model_dir, "libbi")
    
    dir.create(libbi_dir)
  }
  
  ## Make the log file
  logs <- file.path(model_dir, "log.txt")
  con <- file(logs, open = "wt")
  
  ## Send Output to log
  sink(con)
  sink(con, type = "message")
}
  
  # Set the time scale for the model ----------------------------------------


  if (time_scale == "year") {
    time_scale_numeric <- 1
  }else if (time_scale == "month") {
    time_scale_numeric <- 12
  }else{
    stop("Only a yearly (year) or monthly (month) timescale have been implemented.")
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
  
  if (!scale_rate_treat) {
    tb_model_raw <- fix(tb_model_raw, scale_rate_treat = 0)
  }
  
  if (scale_transmission %in% "none") {
    tb_model_raw <- fix(tb_model_raw, beta_df = 1)
  }else if (scale_transmission %in% "young_adult") {
    tb_model_raw <- fix(tb_model_raw, beta_df = 2)
  }
  
  if (scale_nonukborn_mixing %in% "none") {
    tb_model_raw <- fix(tb_model_raw, M_df = 1)
  }else if (scale_nonukborn_mixing %in% "young_adult") {
    tb_model_raw <- fix(tb_model_raw, M_df = 2)
  }
  
  if (!noise) {
    tb_model_raw <- fix(tb_model_raw, noise_switch = 0)
  }
  
  if (!initial_uncertainty) {
    tb_model_raw <- fix(tb_model_raw, initial_uncertainty_switch = 0)
  }
  
  
# Add the proposal block to the model -------------------------------------
 if (!is.null(proposal_param_block)) {
   tb_model_raw <- add_block(tb_model_raw, "proposal_parameter", proposal_param_block)
 }
  
  
  if (!is.null(proposal_initial_block)) {
    tb_model_raw <- add_block(tb_model_raw, "proposal_initial", proposal_initial_block)
  }


  # Allow for interactive debugging -----------------------------------------
  
  if (browse) {
    browser() 
  }


# Set up model input ------------------------------------------------------

## See ?setup_model_input or details 
input <- setup_model_input(run_time, time_scale_numeric)

## Save formated data
if (save_output) {
  saveRDS(input, file.path(data_dir, "input.rds"))
}


# Set up abs data ---------------------------------------------------------

## See ?setup_model_obs for details
obs <- setup_model_obs(years_of_age = years_of_age, age_groups = age_groups,
                       con_age_groups = con_age_groups, spacing_of_historic_tb = spacing_of_historic_tb,
                       years_of_data = years_of_data, aggregated = aggregated_observed)

obs <- obs %>% 
  map(~ filter(., time <= run_time)) %>% 
  map( ~ drop_na(., value))

# Reset obs and input if running with test SIR Model ----------------------
  if (model == "SIR") {
    obs <- NULL
    input <- NULL
  }
  
# Generate data from the model --------------------------------------------
  if (gen_data) {
    
    if (verbose) {
      message("Generating data from model")
    }
    
    tb_data <- bi_generate_dataset(tb_model_raw, end_time = run_time * time_scale_numeric, 
                                   input = input, obs = obs, 
                                   noutputs = floor(run_time / 7),
                                   seed = seed, verbose = libbi_verbose)
    
    if (verbose) {
      message("Summary of generated model data")
      print(tb_data)
    }
    
    
    obs <- bi_read(tb_data, type = "obs")
    
    if (verbose) {
      message("Variables in the generated data")
      print(names(obs))
    }
    
  }
  
  ## Save formated data
  if (save_output) {
    saveRDS(obs, file.path(data_dir, "obs.rds"))
  }
  
  
  # Plot incidence generated by the model -----------------------------------
  if (plot_obs) {
    
    time <- NULL; value <- NULL;
    
    message("Plot number of cases detected each year:")
    if (!is.null(years_of_age) && !is.null(age_groups)) {
      p_cases <- obs$YearlyAgeInc %>% 
        dplyr::filter(time > 0) %>% 
        dplyr::group_by(time) %>% 
        dplyr::summarise(value = sum(value, na.rm = TRUE)) %>% 
        ggplot(aes(x = time, y = value)) +
        geom_point(size = 1.2) +
        geom_line(size = 1.1, alpha = 0.6) +
        theme_minimal() +
        labs(x = "Time",
             y = "Yearly notified cases")
      
      print(p_cases)
      
      if (save_output) {
        ggsave("obs-cases.png", path = plots_dir, dpi = 320)
      }
    }
    
  message("Plot number of cases detected each year, stratified by age group:")
  if (!is.null(years_of_age) && !is.null(age_groups)) {
    p_age_cases <- obs$YearlyAgeInc %>% 
      dplyr::filter(time > 0) %>% 
      ggplot(aes(x = time, y = value)) +
      geom_point(size = 1.2) +
      geom_line(size = 1.1, alpha = 0.6) +
      theme_minimal() +
      labs(x = "Time",
           y = "Yearly notified cases") +
      facet_wrap(~age, scales = "free_y")
    
    print(p_age_cases)
    
    if (save_output) {
      ggsave("obs-age-cases.png", path = plots_dir, dpi = 320)
    }
  }
  }


  # Set the number of particles ---------------------------------------------

  if (is.null(nparticles)) {
  nparticles <- obs %>% 
    map_dbl(nrow) %>%
    sum
  
  nparticles <- 2**ceiling(log(nparticles, 2))
  
  if (verbose) {
    message("Using ", nparticles, " particles based on the number of observed data samples (rounded to the nearest power of 2).")
  }
  
  } 


if (adapt_particles) {
  if (is.null(min_particles)) {
    min_particles <- nparticles
    
    if (verbose) {
      message("Using a minimum of ", min_particles, " particles based on the number of particles specified.")
    }
  }  
  
  if (is.null(max_particles)) {
    max_particles <- min_particles * 4
    
    if (verbose) {
      message("Using a maximum of ", max_particles, " particles based on the minimum number of particles times by 4.")
    }
  }  
}

  # Load the model ----------------------------------------------------------
  
  if (verbose) {
    message("Load the model as a compiled libbi object")
  }

  tb_model <- libbi(tb_model_raw, 
                    input = input, 
                    obs = obs,
                    end_time = run_time * time_scale_numeric, 
                    nparticles = nparticles, 
                    nthreads = nthreads, 
                    debug = libbi_verbose,
                    assert = FALSE,
                    single = TRUE)
  
  if (!is.null(previous_model_path)) {
    message("Replacing the default liBBi model with a previously run model - this may not have the same settings as the current run.")
    tb_model <- read_libbi(previous_model_path)
    
    if (verbose) {
      message("Previous Model: ")
      print(tb_model) 
      print(tb_model$model)
    }
  }
  # Sample from priors ------------------------------------------------------
if (sample_priors) {
  if (verbose) {
    message("Sample priors")
  }
  
  priors <- sample(tb_model, target = "prior", 
                   nsamples = prior_samples)
  
  if (verbose) {
    message("Summary of prior sampling")
    print(priors)
   
  }
  
  if (save_output) {
    if (verbose) {
      message("Save prior samples")
    }
    
    priors %>% 
      bi_read(type = c("param")) %>% 
        map_dfr(~mutate(., distribution = "Prior"), .id = "parameter") %>% 
    saveRDS(file.path(data_dir, "prior-params.rds"))
  }
}  

  

# Optimise the deterministic model ----------------------------------------
if (optim) {
  
  if(verbose) {
    message("Optimising the deterministic model")
  }
  
  tb_model <- tb_model %>% 
    optimise()
  
}
  
# Adapting the number of particles ----------------------------------------
  
  if (adapt_particles) {
    if (verbose) {
      message("Adapting particles starting with ", min_particles, " up to a maximum of ", max_particles)
    }
    
    adapt_mutli_particles <- function(iteration, model, min_particles, max_particles) {
      if (verbose) {
        message("Adapting particles iteration: ", iteration)
        message("Initial sampling using the prior as the proposal.")
      }
      
      model <- model %>% 
        sample(proposal = "prior", nsamples = adapt_part_samples,
               seed = iteration + seed)
      
      if (verbose) {
        message("Starting particle adaption")
      }
      
      model <- adapt_particles(model,
                               min = min_particles, 
                               max = max_particles,
                               quiet = !verbose,
                               target.variance = 5)
      
      particles <- tb_model$options$nparticles
      return(particles)
    }
    
    particles <- map_dbl(1:adapt_part_it, ~  
                           adapt_mutli_particles(., 
                                                 model = tb_model,
                                                 min_particles = min_particles,
                                                 max_particles = max_particles))

    
    if (save_output) {
      saveRDS(particles, file.path(libbi_dir, "particles.rds"))
    }
    
    med_particles <- particles %>% 
      median %>% 
      ceiling
    
    if (verbose) {
      message("Choosing ", med_particles, "based on ", adapt_part_it, " iterations each using ", adapt_part_samples, ".")
      message("The maximum number of particles selected was ", max(particles), " with a minimum of ", min(particles))
      message("The mean number of particles selected was ", mean(particles), " with a standard deviation of ", sd(particles))
    }
    
 
    tb_model$options$nparticles <- med_particles
    nparticles <- med_particles
    
  }
  
  
  # Adapting the proposal ---------------------------------------------------
  
  if (adapt_proposal) {
    if (verbose) {
      message("Adapting proposal with a min acceptance of ", min_acc, " and a maximum acceptance of ", max_acc)
      message("Running for ", adapt_prop_it, " iterations with ", adapt_prop_samples, " samples each time.")
    }
    
    ## Adapt proposal
    if(verbose) {
      message("Taking initial sample to estimate acceptance rate")
    }
    
    tb_model <- tb_model %>% 
      sample(target = "posterior", proposal = "model",
             nsamples = adapt_prop_samples) 
    
    if(verbose) {
      message("Adapting proposal ...")
    }
    
    tb_model <- tb_model %>% 
      adapt_proposal(min = min_acc, max = max_acc,
                     nsamples = adapt_prop_samples,
                     max_iter = adapt_prop_it, adapt = adapt, 
                     scale = adapt_scale, truncate = TRUE, 
                     quiet = !verbose)
    
    acc_rate <- acceptance_rate(tb_model)
    

    
    prop_param_block <- get_block(tb_model$model, "proposal_parameter")
    prop_initial_block <- get_block(tb_model$model, "proposal_initial")
    
    if (verbose) {
      message("Adapted proposal distribution")
    
      print(prop_param_block)
      print(prop_initial_block)
    }
    
    if (save_output) {
      saveRDS(prop_param_block, file.path(libbi_dir, "proposal-param-block.rds"))
      saveRDS(prop_initial_block, file.path(libbi_dir, "proposal-initial-block.rds"))
    }
  }
  

# Set number of moves for SMC2 --------------------------------------------

if (is.null(rejuv_moves)) {
  
  if (!adapt_proposal) {
    if (verbose) {
      message("Taking a 100 samples from the posterior to estimate
              the acceptance rate as proposal adaption has not been run.")
    }
    
    tb_model <- tb_model %>% 
      sample(target = "posterior", proposal = "model",
             nsamples = 1000, 
             verbose = libbi_verbose) 
  }
  
  acc_rate <- acceptance_rate(tb_model)
  
  if (verbose) {
    message("Acceptance rate of ", acc_rate, " after adapting the proposal")
  }

  if(acc_rate < 0.001) {
    if (verbose) {
      message("Acceptance rate is to low (", acc_rate, ") to be tractable. Defaulting to an acceptance rate of 0.001")
    }
    acc_rate <- 0.001
  } 

  
  target_acc <- sample_ess_at
  rejuv_moves <- round(target_acc / acc_rate, digits = 0)
  rejuv_moves <- ifelse(rejuv_moves < 1, 1, rejuv_moves)
  
  if (verbose) {
    message("Using ", rejuv_moves, " rejuvernation moves in order to target at least a ", round(100*target_acc, 0), "% acceptence rate for each rejuvernation sample.")
  }
}  
  
  
  # Model fitting ------------------------------------------------------
  
  
  if (fit) {
   
    if (verbose) {
      if (time_for_resampling != 0) { 
        message("As the time for resampling has been specified rejuvernation will happen after each SMC step and will take as long as has been allocated regardless of the acceptance rate.")
        }
      message("Fitting using SMC-SMC")
    }
    
    tb_model <- tb_model %>% 
      sample(target = "posterior",
             proposal = "model",
             nparticles = nparticles,
             nsamples = posterior_samples,
             sampler = "sir", 
             `sample-ess-rel` = sample_ess_at,
             nmoves = rejuv_moves,
             tmoves = time_for_resampling * 60,
             thin = thin,
             verbose = fitting_verbose)
    
    
    if (verbose) {
      message("Summary of fitted model")
      print(tb_model)
    }
    

# Evaluate run ------------------------------------------------------------
    
    burn <- floor(posterior_samples / thin * burn_prop)
    
    dic <- DIC(tb_model, burn = burn)
    
    message("Model DIC: ", dic)
    
    if (save_output) {
      saveRDS(dic, file = file.path(libbi_dir, "dic.rds"))
    }

  }
  

# Predict states ----------------------------------------------------------
  if (!fit && !adapt_particles && !adapt_proposal && sample_priors) {
    tb_model <- priors
  }
  
if (pred_states ) {
  
  if (verbose) {
    message("Predicting states based on posterior sample up to 2040")
  }
  
  ## Predicting states for all times from intialisation to end time (for the standard model in 2040).
  tb_model <- sample_obs(tb_model, 
                         end_time = (36 + run_time) * time_scale_numeric, 
                         noutputs = (36 + run_time) * time_scale_numeric)
  
}  
  

# Save model --------------------------------------------------------------
  
  if (save_output) {
    save_libbi(tb_model, file.path(libbi_dir, "posterior"))
  }

# Model report ------------------------------------------------------------

  
  if (save_output && reports && fit) {
    
    if (verbose) {
      ## Assumes the function is being run at the root of the project (may need refinement)
      message("Creating model report")
    }
    
    ## Model report
    model_report <- "./vignettes/model_report.Rmd"
    future_scenarios <- "./vignettes/future_scenarios.Rmd"
    
    report_dir <- file.path(model_dir, "reports")
    dir.create(report_dir)
    
    rmarkdown::render(model_report, output_format = "html_document",
                      output_dir = report_dir, 
                      knit_root_dir = "..",
                      params =  list(model_dir =  model_dir))
    
    
    rmarkdown::render(future_scenarios, output_format = "html_document", 
                      output_dir = report_dir, 
                      knit_root_dir = "..",
                      params =  list(model_dir =  model_dir))
    
  }
  
  if (save_output) {
    sink(type = "message")
    sink() 
  }
  
  if (!exists("tb_model")) {
    tb_model <- NULL
  }
  return(tb_model)
}
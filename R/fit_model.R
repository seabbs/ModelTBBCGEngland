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
#' is used as the initial particle number.
#' @param nthreads Numeric, defaults to 4. The number of parallel jobs to run. The most efficient option is likely to be to match the 
#' number of cores available.
#' @param adapt_particles Logical, defaults to \code{FALSE}. Should the number of particles be adapted.
#' @param adapt_part_samples Numeric, the number of samples to take from the priors when adapting the number of particles. Defaults to 1000.
#' @param adapt_part_it Numeric, the number of iterations to use for adapting the number of particles. Defaults to 10.
#' @param min_particles Numeric, defaults to 4. The starting number of particles to use for adpation. If \code{NULL} uses \code{nparticles}.
#' @param max_particles Numeric, defaults to 16. The maximum number of particles to use for apation. If \code{NULL} then defaults to \code{4 * min_particles}
#' @param proposal_param_block A character string containing a LiBBi proposal parameter block. Defaults to \code{NULL} in 
#' which case the LiBBi defaults will be used.
#' @param proposal_initial_block A character string containing a LiBBi proposal initial block. Defaults to \code{NULL} in 
#' which case the LiBBi defaults will be used.
#' @param adapt_proposal  Logical, defaults to \code{TRUE}. Should the proposal be adjusted to improve the acceptance rate.
#' @param adapt_part_samples Numeric, the number of samples to take from the priors when adapting the proposal distribution. Defaults to 1000.
#' @param adapt_part_it Numeric, defaults to 10. The number of iterations to use for adapting the proposal.
#' @param adapt Character string, defaults to "both". The type of adaption to use for the proposal see \code{rbi.helpers::adapt_proposal} for details.
#' @param adapt_scale Numeric, defaults to 2. The scale to use to adapt the proposal.
#' @param min_acc Numeric, defaults to 0.05. The minimum target acceptance rate.
#' @param max_acc Numeric, defaults to 0.4. The maximum target acceptance rate.
#' @param fit Logical, defaults to \code{TRUE}. Should the model be fitted with 1000 samples extracted.
#' @param thin Numeric, the thinning interval to use between posterior samples. May reduce the correlation between samples and reduces memory.
#' @param burn_prop Numeric (between 0 and 1). The proportion of the chain to discard as burn-in.
#' @param sample_ess_at Numeric defaults to 0.8. The thresold of the effective sample size (ess) at which to rejuvernate the particles.
#' @param rejuv_moves Numeric, defaults to \code{NULL}. The number of MCMC samples to take for each rejuvernation step. If \code{NULL} is set so that
#' the acceptance rate of each rejuvernation step is at least 20%. This is based on the estimated acceptance rate after adaption. If proposal adaption is not used
#' then this defaults to 1.
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
#' @param pred_states Logical defaults to \code{TRUE}. Should states be predicted over all time (from model initialisation to 35 years ahead of final run time). 
#' If set to \code{FALSE} states will only be estimated for times with observed data points.
#' @param seed Numeric, the seed to use for random number generation.
#' @param reports Logical, defaults to \code{TRUE}. Should model reports be generated. Only enabled when \code{save_output = TRUE}.
#' @return A LibBi model object based on the inputed test model.
#' @export
#' @inheritParams setup_model_obs
#' @importFrom rbi fix bi_model sample bi_read bi_generate_dataset libbi get_block
#' @import rbi.helpers 
#' @import ggplot2
#' @importFrom dplyr filter mutate select vars arrange count rename
#' @importFrom tidyr drop_na
#' @importFrom stats runif time
#' @importFrom utils str
#' @importFrom graphics plot
#' @importFrom purrr map map_dbl
#' @importFrom tibble tibble
#' @importFrom stringr str_replace
#' @importFrom rmarkdown render
#' @examples
#' 
#' ## Function code:
#' fit_model
fit_model <- function(model = "BaseLineModel", previous_model_path = NULL, gen_data = TRUE, 
                      run_time = 73, time_scale = "year", plot_obs = TRUE,
                      sample_priors = TRUE, prior_samples = 1000, nparticles = NULL,
                      adapt_particles = TRUE, adapt_part_samples = 1000, adapt_part_it = 10, 
                      min_particles = NULL, max_particles = NULL,
                      proposal_param_block = NULL, proposal_initial_block = NULL, 
                      adapt_proposal = TRUE, adapt_prop_samples = 100, adapt_prop_it = 3, adapt = "both",
                      adapt_scale = 2, min_acc = 0.04, max_acc = 0.4,
                      fit = FALSE, posterior_samples = 10000, thin = 1, burn_prop = 0, sample_ess_at = 0.8,
                      rejuv_moves = NULL, nthreads = parallel::detectCores(), verbose = TRUE, libbi_verbose = FALSE, 
                      fitting_verbose = TRUE, pred_states = TRUE, browse = FALSE,
                      const_pop = FALSE, no_age = FALSE, no_disease = FALSE, scale_rate_treat = TRUE, years_of_age = c(2000, 2004),
                      spacing_of_historic_tb = 5, noise = TRUE,
                      save_output = FALSE, dir_path = NULL, dir_name = NULL, reports = TRUE,
                      seed = 1234) {


# Util functions ----------------------------------------------------------

  ##Plot and save
  plot_obj <- function(obj, p_param, append_name = NULL, save = FALSE) {
    if (!is.null(p_param[[obj]]))
    {
      message("Plotting: ", obj)
      
      plot <- p_param[[obj]] + theme_minimal()
      
      print(plot)
      
      if (!is.null(append_name)) {
        obj <- paste0(append_name, "-", obj)
      }
      
      if (save) {
        ggsave(paste0(obj, ".png"), plot, path = plots_dir, dpi = 320)
      }

    }
  }

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
  con <- file(logs)
  
  ## Send Output to log
  sink(con, append = TRUE)
  sink(con, type = "message", append = TRUE)
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
  
  if (!noise) {
    tb_model_raw <- fix(tb_model_raw, noise_switch = 0)
    nparticles < - 1
    adapt_particles <- FALSE
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
obs <- setup_model_obs(years_of_age = years_of_age, spacing_of_historic_tb = spacing_of_historic_tb)

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
    
    if (save_output) {
      ggsave("obs-hist-pul-cases.png", path = plots_dir, dpi = 320)
    }
    
    message("Plot number of cases detected each year:")
    p_cases <- obs$YearlyInc %>% 
      dplyr::filter(time > 0) %>% 
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
    
  message("Plot number of cases detected each year, stratified by age group:")
  if (!is.null(years_of_age)) {
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
    map(~ filter(., time <= run_time)) %>% 
    map( ~ drop_na(., value)) %>% 
    map_dbl(nrow) %>%
    sum
  
  if (verbose) {
    message("Using ", nparticles, " particles based on the number of observed data samples.")
  }
  
  } 

if (is.null(min_particles)) {
  min_particles <- round(nparticles / 2, digits = 0)
  
  if (verbose) {
    message("Using a minimum of ", min_particles, " particles based on the number of observed data samples divided by two.")
  }
}  
  
if (is.null(max_particles)) {
  max_particles <- round(min_particles * 4, digits = 0)
  
  if (verbose) {
    message("Using a maximum of ", max_particles, " particles based on the minimum number of particles times by 4.")
  }
}  
  # Load the model ----------------------------------------------------------
  
  if (verbose) {
    message("Load the model as a compiled libbi object")
  }

  tb_model <- libbi(tb_model_raw, 
                    input = input, 
                    obs = obs,
                    noutputs = run_time,
                    end_time = run_time * time_scale_numeric, 
                    nparticles = nparticles, nthreads = nthreads, 
                    verbose = libbi_verbose,
                    seed = seed)
  
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
  
  priors <- sample(tb_model, target = "prior", nsamples = prior_samples, sample_obs = TRUE)
  
  if (verbose) {
    message("Summary of prior sampling")
    print(priors)
    message("Load priors")
  }
  
  
  if (verbose) {
    message("Plot priors")
  
    p_prior <- plot(priors , plot = FALSE)
    
    objects <- c("states", "densities", "traces", "correlations", "noises", "likelihoods")
    
    map(objects, ~ plot_obj(., p_prior, append_name = "prior", save = save_output))
      
    if (save_output) {
      saveRDS(p_prior$data, file.path(data_dir, "prior-params.rds"))
    }

  }
}  

  # Adapting the number of particles ----------------------------------------
  
  if (adapt_particles) {
    if (verbose) {
      message("Adapting particles starting with ", min_particles, " up to a maximum of ", max_particles)
    }
    
    ## Add all variables as outputs to the model.
    lim_out_model <- tb_model$model
    tb_model$model <- everything_from_model(tb_model$model)

    adapt_mutli_particles <- function(iteration, libbi) {
      if (verbose) {
        message("Adapting particles iteration: ", iteration)
        message("Initial sampling using the prior as the proposal.")
      }
      
      tb_model <- tb_model %>% 
        sample(proposal = "prior", nsamples = adapt_part_samples, verbose = fitting_verbose,
               options = list(with="transform-initial-to-param"), seed = iteration + seed)
      
      if (verbose) {
        message("Starting particle adaption")
      }
      
      tb_model <- tb_model %>% 
        adapt_particles(min = min_particles, max = max_particles, 
                        quiet = !verbose,
                        verbose = libbi_verbose)
      
      particles <- tb_model$options$nparticles
      return(particles)
    }
    
    particles <- map_dbl(1:adapt_part_it, ~  adapt_mutli_particles(., tb_model))

    
    if (save_output) {
      saveRDS(particles, file.path(libbi_dir, "particles.rds"))
    }
    
    med_particles <- particles %>% 
      median %>% 
      ceiling
    
    if (verbose) {
      message("Choosing ", med_particles, "based on ", adapt_part_it, " iterations each using ", adapt_part_samples, ".")
      message("The maximum number of particles selected was ", max(particles), " with a minimum of", min(particles))
      message("The mean number of particles selected was ", mean(particles), " with a standard deviation of ", sd(particles))
    }
    
 
    tb_model$model <- lim_out_model
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
      sample(proposal = "model", nsamples = adapt_prop_samples, verbose = libbi_verbose) 
    
    if(verbose) {
      message("Adapting proposal ...")
    }
    
    tb_model <- tb_model %>% 
      adapt_proposal(min = min_acc, max = max_acc,
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
  
  if (adapt_proposal) {
    acc_rate <- acceptance_rate(tb_model)
    if (verbose) {
      message("Acceptance rate of ", acc_rate, " after adapting the proposal")
    }
  }else{
    acc_rate <- 0.02
    if (verbose) {
      message("Acceptance rate of ", acc_rate, " assumed by default as proposal not adapted.")
    }
  }

  if(acc_rate < 0.0002) {
    if (verbose) {
      message("Acceptance rate is to low (", acc_rate, ") to be tractable. Defaulting to an acceptance rate of 0.0002 (leading to 1000 moves per particle).")
    }
    acc_rate <- 0.002
  } 

  
  target_acc <- 0.2
  rejuv_moves <- round(target_acc / acc_rate, digits = 0)
  rejuv_moves <- ifelse(rejuv_moves < 1, 1, rejuv_moves)
  
  if (verbose) {
    message("Using ", rejuv_moves, " rejuvernation moves in order to target at least a 20% acceptence rate of the MCMC sampler.")
  }
}  
  
  
  # Model fitting ------------------------------------------------------
  
  
  if (fit) {
   
    if (verbose) {
      message("Fitting using PMCMC")
    }
    
    tb_model <- tb_model %>% 
      sample(target = "posterior",
             proposal = "model", sample_obs = TRUE, 
             nparticles = nparticles,
             nsamples = posterior_samples,
             options = list("sampler" = "sir", 
                            "adapter" = "global",
                            "sample-ess-rel" = sample_ess_at,
                            "nmoves" = rejuv_moves),
             thin = thin,
             verbose = fitting_verbose)
    
    
    if (verbose) {
      message("Summary of fitted model")
      print(tb_model)
    }
    
# Plot posterior ----------------------------------------------------------

    if (verbose && fit) {
      message("Plot posterior")
      
      p_posterior <- plot(tb_model, plot = FALSE)
      
      objects <- c("states", "densities", "traces", "correlations", "noises", "likelihoods")
      
      map(objects, ~ plot_obj(., p_posterior, append_name = "posterior", save = save_output))
      
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
  
if (pred_states) {
  
  if (verbose) {
    message("Predicting states based on posterior sample up to 2040")
  }
  
  ## Predicting states for all times from intialisation to end time (for the standard model in 2040).
  tb_model <- predict(tb_model, end_time = (36 + run_time) * time_scale_numeric, 
                      noutputs = (36 + run_time) * time_scale_numeric,
                      verbose = FALSE,
                      sample_obs = TRUE)
  
}  
  

# Save model --------------------------------------------------------------
  
  if (save_output) {
    save_libbi(tb_model, file.path(libbi_dir, "posterior"))
  }

# Model report ------------------------------------------------------------

  
  if (save_output && reports) {
    
    if (verbose) {
      ## Assumes the function is being run at the root of the project (may need refinement)
      message("Creating model report")
    }
    
    rm(tb_model)
    
    report <- "./vignettes/model_report.Rmd"
    report_dir <- file.path(model_dir, "reports")
    dir.create(report_dir)
    
    rmarkdown::render(report, output_format = "html_document", output_dir = report_dir, 
                      knit_root_dir = "..", params =  list(model_dir =  model_dir))
    
    
    
    
    sink(file = NULL) 
    sink(file = NULL, type = "message")
  }
  
  if (!exists("tb_model")) {
    tb_model <- NULL
  }
  return(tb_model)
}
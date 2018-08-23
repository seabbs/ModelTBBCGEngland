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
#' @param nparticles Numeric, the initial number of particles to use in the particle filters. Defaults to \code{NULL}, in which case the number of outputs 
#' is used as the initial particle number.
#' @param adaption_samples s Numeric, the number of samples to take from the priors when adapting the number of particles and the proposal distribution. Defaults to 1000.
#' @param nthreads Numeric, defaults to 4. The number of parallel jobs to run. The most efficient option is likely to be to match the 
#' number of cores available.
#' @param adapt_particles Logical, defaults to \code{FALSE}. Should the number of particles be adapted.
#' @param min_particles Numeric, defaults to 4. The starting number of particles to use for adpation.
#' @param max_particles Numeric, defaults to 16. The maximum number of particles to use for apation.
#' @param proposal_param_block A character string containing a LiBBi proposal parameter block. Defaults to \code{NULL} in 
#' which case the LiBBi defaults will be used.
#' @param proposal_initial_block A character string containing a LiBBi proposal initial block. Defaults to \code{NULL} in 
#' which case the LiBBi defaults will be used.
#' @param adapt_proposal  Logical, defaults to \code{TRUE}. Should the proposal be adjusted to improve the acceptance rate.
#' @param adapt_scale Numeric, defaults to 2. The scale to use to adapt the proposal.
#' @param min_acc Numeric, defaults to 0.05. The minimum target acceptance rate.
#' @param max_acc Numeric, defaults to 0.4. The maximum target acceptance rate.
#' @param fit Logical, defaults to \code{TRUE}. Should the model be fitted with 1000 samples extracted.
#' @param thin Numeric, the thinning interval to use between posterior samples. May reduce the correlation between samples and reduces memory.
#' @param burn_prop Numeric (between 0 and 1). The proportion of the chain to discard as burn-in.
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
#' @importFrom rbi fix bi_model sample bi_read bi_generate_dataset libbi get_block save_libbi read_libbi
#' @import rbi.helpers 
#' @import ggplot2
#' @importFrom dplyr filter mutate select vars arrange count rename
#' @importFrom stats runif time
#' @importFrom utils str
#' @importFrom graphics plot
#' @importFrom purrr map
#' @importFrom tibble tibble
#' @importFrom tidyr unnest
#' @importFrom stringr str_replace
#' @examples
#' 
#' ## Function code:
#' fit_model
fit_model <- function(model= "BaseLineModel", previous_model_path = NULL, gen_data = TRUE, run_time = 73, time_scale = "year", plot_obs = TRUE,
                      sample_priors = TRUE, prior_samples = 1000, nparticles = NULL, adaption_samples = 1000, adapt_particles = FALSE,
                      min_particles = 4, max_particles = 16, proposal_param_block = NULL, proposal_initial_block = NULL,
                      adapt_proposal = FALSE, adapt_scale = 2, min_acc = 0.05, max_acc = 0.4, 
                      fit = FALSE, posterior_samples = 10000, thin = 0, burn_prop = 0, 
                      nthreads = parallel::detectCores(), verbose = TRUE, libbi_verbose = FALSE, browse = FALSE,
                      const_pop = FALSE, no_age = FALSE, no_disease = FALSE,
                      save_output = FALSE, dir_path = NULL, dir_name = NULL) {


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
  file(logs)
  
  ## Send Output to log
  sink(logs, type=c("output", "message"), append = TRUE, split = TRUE)
}
  
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
  
  ## Non UK born cases in 2000 - used to estimate historic non UK born cases
  NUKCases2000 <- nonukborn_p_cases %>% 
    filter(time == time_scale_numeric * (2005 - 1931)) %>% 
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

## Save formated data
saveRDS(input, file.path(data_dir, "input.rds"))

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
                                   noutputs =  10,
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
  
  ## Save formated data
  saveRDS(obs, file.path(data_dir, "obs.rds"))
  
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

  
  # Load the model ----------------------------------------------------------
  
  if (verbose) {
    message("Load the model as a compiled libbi object")
  }

  tb_model <- libbi(tb_model_raw, 
                    input = input, 
                    obs = obs,
                    end_time = run_time * time_scale_numeric, 
                    noutputs = run_time * time_scale_numeric,
                    nparticles = nparticles, nthreads = nthreads, 
                    verbose = libbi_verbose)
  
  if (!is.null(previous_model_path)) {
    message("Replacing the default liBBi model with a previously run model - this may not have the same settings as the current run.")
    tb_model <- read_libbi(previous_model_path)
    
    if (verbose) {
      message("Previous Model: ")
      print(tb_model) 
    }
  }
  # Sample from priors ------------------------------------------------------
if (sample_priors) {
  if (verbose) {
    message("Sample priors")
  }
  
  tb_model <- sample(tb_model, target="prior", nsamples = prior_samples)
  
  if (verbose) {
    message("Summary of prior sampling")
    print(tb_model )
    message("Load priors")
  }
  
  
  if (verbose) {
    message("Plot priors")
  
    p_prior <- plot(tb_model , plot = FALSE)
    
    objects <- c("states", "densities", "traces", "correlations", "noises", "likelihoods")
    
    map(objects, ~ plot_obj(., p_prior, append_name = "prior", save = save_output))
      
    if (save_output) {
      saveRDS(p_prior$data, file.path(data_dir, "prior-params.rds"))
    }

  }
  
  if (save_output) {
    save_libbi(tb_model, file.path(libbi_dir, "priors.rds"))
  }
  
  
}  

  # Adapting the number of particles ----------------------------------------
  
  if (adapt_particles) {
    if (verbose) {
      message("Adapting particles starting with ", min_particles, " up to a maximum of ", max_particles)
    }
    
    tb_model <- adapt_particles(tb_model,
                                nsamples = adaption_samples, 
                                min = min_particles, max = max_particles, 
                                quiet = !verbose,
                                verbose = libbi_verbose,
                                target = "prior")
    
    particles <- tb_model$options$nparticles
    
    if (verbose) {
      message(particles, " has been selected as meeting the variance criteria")
    }
    
    if (save_output) {
      saveRDS(particles, file.path(libbi_dir, "particles.rds"))
    }
    
  }
  
  
  # Adapting the proposal ---------------------------------------------------
  
  if (adapt_proposal) {
    if (verbose) {
      message("Adapting proposal with a min acceptance of ", min_acc, " and a maximum acceptance of ", max_acc)
    }
    
    tb_model <- adapt_proposal(tb_model, min = min_acc, max = max_acc,
                               max_iter = 5, nsamples = adaption_samples,
                               adapt = c("both"), scale = adapt_scale,
                               truncate = TRUE, quiet = !verbose,
                               verbose = libbi_verbose,
                               target= "prior")
    
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
  
  
  # Model fitting ------------------------------------------------------
  
  
  if (fit) {
    if (verbose) {
      message("Fitting using PMCMC")
    }
    
    
    tb_model <- tb_model %>% 
      sample(target="posterior", sample_obs = TRUE, 
             nsamples = posterior_samples,
             thin = thin,
             verbose = libbi_verbose)
    
    
    if (verbose) {
      message("Summary of fitted model")
      print(tb_model)
    }
    
    if (save_output) {
      save_libbi(tb_model, file.path(libbi_dir, "posterior.rds"))
    }
    

# Plot posterior ----------------------------------------------------------

    if (verbose) {
      message("Plot posterior")
      
      p_posterior <- plot(tb_model, plot = FALSE)
      
      objects <- c("states", "densities", "traces", "correlations", "noises", "likelihoods")
      
      map(objects, ~ plot_obj(., p_posterior, append_name = "posterior", save = save_output))
      
      if (save_output) {
        saveRDS(p_posterior$data, file.path(data_dir, "posterior-params.rds"))
      }
      
    }
    
    if (save_output) {
      save_libbi(tb_model, file.path(libbi_dir, "posterior.rds"))
    }
    

# Evaluate run ------------------------------------------------------------
    
    burn <- floor(posterior_samples / thin * burn_prop)
    
    dic <- DIC(tb_model, burn = burn)
    
    message("Model DIC: ", dic)
    
    if (save_output) {
      saveRDS(dic, file = file.path(libbi_dir, "dic.rds"))
    }

  }
  

  if (save_output) {
    sink(file = NULL) 
  }
  
  return(tb_model)
}
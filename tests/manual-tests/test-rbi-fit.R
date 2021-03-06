library('rbi')
library('rbi.helpers')
library('ModelTBBCGEngland')

## Should particles be adapted
gen_data <- FALSE
sample_priors <- TRUE
adapt_part <- FALSE
adapt_prop <- FALSE
sample_post <- FALSE
use_sir_sampling <- TRUE 
pred_sample <- TRUE
verbose <- TRUE
save_results <- TRUE
det_optim <- FALSE
model <- "BaseLineModel" ##"BaseLineModel"
noise <- FALSE
initial_uncertainty <- FALSE
trans_prob_freedom <- "young_adult"
free_non_uk_mixing <- TRUE
show_conf <- FALSE

if (use_sir_sampling) {
  sample_post <- FALSE
}

## Need to preload input
input <- setup_model_input(run_time = 73, time_scale_numeric = 1)
obs <- setup_model_obs(years_of_age = c(2000:2004), 
                       age_groups = 0:11,
                       aggregated = FALSE,
                       years_of_data = c(2000:2004))

# Load model file ---------------------------------------------------------

model_file <- system.file(package="ModelTBBCGEngland", paste0("bi/", model, ".bi")) # get full file name from package
tb_model_raw <- bi_model(model_file) # load model

if (model == "BaseLineModel") {
  tb_model_raw <- tb_model_raw %>% 
    fix(no_disease = 0, timestep = 1)
}

if (!noise) {
  tb_model_raw <- tb_model_raw %>% 
    fix(noise_switch = 0)
}

if (trans_prob_freedom %in% "none") {
  tb_model_raw <- fix(tb_model_raw, beta_df = 1)
}else if (trans_prob_freedom %in% "young_adult") {
  tb_model_raw <- fix(tb_model_raw, beta_df = 2)
}


if (free_non_uk_mixing) {
  tb_model_raw <- fix(tb_model_raw, M_df = 2)
}

if (!initial_uncertainty) {
  tb_model_raw <- fix(tb_model_raw, initial_uncertainty_switch = 0)
}


# Generate a simulated dataset --------------------------------------------

if (gen_data) {

  tb_data <- bi_generate_dataset(tb_model_raw, end_time=73, noutputs=4, seed=12345678, input = input)
  
  
  obs <- tb_data
}


# Set up Libbi model ------------------------------------------------------

model <- libbi(tb_model_raw, 
              nsamples = 1000, end_time = 73,
              nparticles = 4, obs = obs, 
              input = input, seed=1234,
              nthreads = 15,
              single = TRUE,
              assert = FALSE)

# Sample priors -----------------------------------------------------------

if (sample_priors) {
  prior <- sample(model, target = "prior", verbose = TRUE,
                  nsamples = 1000, noutputs = 73)
  
  prior <- predict(prior, with = "transform-obs-to-state") 
}


if (sample_priors) {
 model_plots(model = prior, uncertainty = TRUE)
}

# Optimise deterministic model --------------------------------------------

if (det_optim) {
  model <- model %>% 
    optimise(verbose = TRUE)
}

# Adapt particles ---------------------------------------------------------

if (adapt_part) {

adapted <- adapt_particles(model, min = 4, 
                           max = 64, 
                           nsamples = 1000,
                           target.variance = 1,
                           verbose = FALSE)

adapted$options$nparticles
}else{
  adapted <- model
}


# Adapt proposal ----------------------------------------------------------

if (adapt_prop) {
  
  adapted <- adapt_proposal(adapted, min = 0.05,
                            max = 0.15, 
                            adapt = "shape",
                            max_iter = 5,
                            nsamples = 1000, 
                            verbose = TRUE)
  
  get_block(adapted$model, "proposal_parameter")
}


# Sample posterior using PMCMC --------------------------------------------

if (sample_post) {

  posterior <- rbi::sample(adapted,
                      target = "posterior",
                      proposal = "model",
                      nsamples = 1000,
                      thin = 1, verbose = TRUE)

}else{
  posterior <- adapted
}



# Sample posterior using SIR/SMC2 -----------------------------------------


if (use_sir_sampling) {
  posterior_smc <- sample(posterior, target = "posterior", 
                          nsamples = 10000, 
                          sampler = "sir", 
                          nmoves = 5, 
                          `sample-ess-rel` = 0.1,
                          thin = 1,
                          verbose = TRUE)
  
  plot_param(posterior_smc, prior_params = prior)
  
}else{
  posterior_smc <- posterior
}     

# Predict states for all times. -------------------------------------------

if (pred_sample) {
  posterior_smc <- predict(posterior_smc, debug = FALSE,
                           with = "transform-obs-to-state",
                           noutputs = 73)
  
  model_plots(posterior_smc, prior = prior, uncertainty = TRUE)

}

if (save_results) {
ModelTBBCGEngland::save_libbi(posterior_smc, "rbi-test-model")
}
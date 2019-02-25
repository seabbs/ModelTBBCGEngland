library('rbi')
library('rbi.helpers')
library('ModelTBBCGEngland')

## Should particles be adapted
gen_data <- FALSE
sample_priors <- TRUE
adapt_part <- TRUE
adapt_prop <- FALSE
sample_post <- TRUE
use_sir_sampling <- TRUE
pred_sample <- TRUE
verbose <- TRUE
save_results <- FALSE
det_optim <- TRUE
model <- "BaseLineModel" ##"BaseLineModel"

if (use_sir_sampling) {
  sample_post <- FALSE
}

## Need to preload input
input <- setup_model_input(run_time = 73, time_scale_numeric = 1)
obs <- setup_model_obs(years_of_data = 2000,
                      years_of_age = NULL, 
                      con_age_groups = c("children", "older adults"),
                      spacing_of_historic_tb = 10)

# Load model file ---------------------------------------------------------


model_file <- system.file(package="ModelTBBCGEngland", paste0("bi/", model, ".bi")) # get full file name from package
SIRmodel <- bi_model(model_file) # load model

if (model == "BaselineModel") {
  SIRmodel <- fix(SIRmodel, noise = 0) %>% 
    fix(scale_rate_treat = 0)
}

# Generate a simulated dataset --------------------------------------------

if (gen_data) {

  SIRdata <- bi_generate_dataset(SIRmodel, end_time=73, noutputs=4, seed=12345678, input = input)
  
  
  obs <- SIRdata
}


# Set up Libbi model ------------------------------------------------------

model <- libbi(SIRmodel, 
              nsamples = 1000, end_time = 73,
              nparticles = 4, obs = obs, 
              input = input, seed=1234,
              nthreads = 4,
              single = TRUE,
              assert = FALSE)


# Optimise deterministic model --------------------------------------------

if (det_optim) {
  model <- model %>% 
    optimise()
}
# Sample priors -----------------------------------------------------------

if (sample_priors) {
  prior <- sample(model, target = "prior")
}

# Run mcmc using the prior as the proposal --------------------------------

if (adapt_part || adapt_prop) {
  
  bi_prior <- sample(model, proposal="prior", nsamples = 100, verbose = TRUE)
}else{
  bi_prior <- model
}

# Adapt particles ---------------------------------------------------------

if (adapt_part) {
  
tmp_model <- bi_prior$model
bi_prior$model <- everything_from_model(tmp_model)

adapted <- adapt_particles(bi_prior, min = 1, max = 16)

bi_prior$model <- tmp_model

adapted$options$nparticles
}else{
  adapted <- bi_prior
}


# Adapt proposal ----------------------------------------------------------

if (adapt_prop) {
  adapted <- adapt_proposal(adapted, min=0.1, max=0.4, adapt = "size", 
                            scale = 2, max_iter = 5, nsamples = 100, verbose = TRUE)
  
  get_block(adapted$model, "proposal_parameter")
}

if (save_results) {
  ModelTBBCGEngland::save_libbi(adapted, "rbi-test-adapted")
}

# Sample posterior using PMCMC --------------------------------------------

if (sample_post) {
  tic()
  posterior <- rbi::sample(adapted,
                      target = "posterior",
                      proposal = "model",
                      nsamples = 1000,
                      thin = 1, verbose = TRUE)
  toc()
  p <- plot(posterior, plot = FALSE)
  p 
}else{
  posterior <- adapted
}


# Sample posterior using SIR/SMC2 -----------------------------------------


if (use_sir_sampling) {
  posterior <- sample(posterior, target = "posterior", 
                      nsamples = 100, 
                      sampler = "sir", 
                      adapter = "global",
                      tmoves =  0,
                      nmoves = 10,
                      `sample-ess-rel` = 0.1,
                      thin = 1,
                      verbose = TRUE)
  
  p <- plot(posterior, plot = FALSE)
  
  p
}else{
  posterior <- adapted
}


# Predict states for all times. -------------------------------------------

if (pred_sample) {
  posterior <- predict(posterior, end_time = 120, noutputs = 120, debug = FALSE)
  
  p <- plot(posterior, plot = FALSE)
  
  p
}

if (save_results) {
ModelTBBCGEngland::save_libbi(posterior, "rbi-test-model")
}
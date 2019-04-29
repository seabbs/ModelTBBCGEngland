library('rbi')
library('rbi.helpers')
library('ModelTBBCGEngland')

## Should particles be adapted
gen_data <- FALSE
sample_priors <- TRUE
adapt_part <- FALSE
adapt_prop <- FALSE
sample_post <- FALSE
use_sir_sampling <- FALSE
pred_sample <- FALSE
verbose <- TRUE
save_results <- TRUE
det_optim <- TRUE
model <- "BaseLineModel" ##"BaseLineModel"

if (use_sir_sampling) {
  sample_post <- FALSE
}

## Need to preload input
input <- setup_model_input(run_time = 73, time_scale_numeric = 1)
obs <- setup_model_obs(years_of_data = 2000:2004,
                       years_of_age = c(2000:2004), 
                       con_age_groups = c("children", "older adults"),
                       spacing_of_historic_tb = 10)

# Load model file ---------------------------------------------------------


model_file <- system.file(package="ModelTBBCGEngland", paste0("bi/", model, ".bi")) # get full file name from package
SIRmodel <- bi_model(model_file) # load model

if (model == "BaseLineModel") {
  SIRmodel <- SIRmodel %>% 
    fix(no_disease = 0, timestep = 1)
}

# Generate a simulated dataset --------------------------------------------

if (gen_data) {

  SIRdata <- bi_generate_dataset(SIRmodel, end_time=73, noutputs=4, seed=12345678, input = input)
  
  
  obs <- SIRdata
}


# Set up Libbi model ------------------------------------------------------

model <- libbi(SIRmodel, 
              nsamples = 1000, end_time = 73,
              nparticles = 256, obs = obs, 
              input = input, seed=1234,
              nthreads = 16,
              single = TRUE,
              assert = FALSE)

# Sample priors -----------------------------------------------------------

if (sample_priors) {
  prior <- sample(model, target = "prior", verbose = TRUE, nsamples = 1000, noutputs = 73)
}

# Optimise deterministic model --------------------------------------------

if (det_optim) {
  model <- model %>% 
    optimise(verbose = TRUE)
}

# Adapt particles ---------------------------------------------------------

if (adapt_part) {
adapted <- rbi::sample(model,
                       target = "posterior",
                       proposal = "model",
                       nsamples = 250,
                       verbose = TRUE)

adapted <- adapt_particles(adapted, min = 256, max = 1024, nsamples = 100,
                           target.variance = 5)

adapted$options$nparticles
}else{
  adapted <- model
}


# Adapt proposal ----------------------------------------------------------

if (adapt_prop) {
  
  adapted <- adapt_proposal(adapted, min = 0.1,
                            max = 0.2, 
                            adapt = "size",
                            max_iter = 2,
                            nsamples = 100, 
                            verbose = TRUE)
  
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
                      nsamples = 250,
                      thin = 1, verbose = TRUE)
  toc()
}else{
  posterior <- adapted
}


# Sample posterior using SIR/SMC2 -----------------------------------------


if (use_sir_sampling) {
  posterior_smc <- sample(posterior, target = "posterior", 
                      nsamples = 1000, 
                      sampler = "sir", 
                      nmoves = 1,
                      `sample-ess-rel` = 0.1 ,
                      thin = 1,
                      verbose = TRUE)

}else{
  posterior <- adapted
}


# Predict states for all times. -------------------------------------------

if (pred_sample) {
  posterior_smc <- predict(posterior_smc, end_time = 120, noutputs = 120, debug = FALSE)
}

if (save_results) {
ModelTBBCGEngland::save_libbi(posterior_smc, "rbi-test-model")
}
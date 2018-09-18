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
verbose <- FALSE
save_results <- FALSE
model <- "BaseLineModel"

if (use_sir_sampling) {
  sample_post <- FALSE
}

## Need to preload input
input <- setup_model_input(run_time = 73, time_scale_numeric = 1)
obs <- setup_model_obs(years_of_age = c(2000, 2004))

# Load model file ---------------------------------------------------------


model_file <- system.file(package="ModelTBBCGEngland", paste0("bi/", model, ".bi")) # get full file name from package
SIRmodel <- bi_model(model_file) # load model

if (model == "BaselineModel") {
  SIRmodel <- fix(SIRmodel, noise = 0) %>% 
    fix(scale_rate_treat = 0)
}

# Generate a simulated dataset --------------------------------------------

if (gen_data) {

  SIRdata <- bi_generate_dataset(SIRmodel, end_time=73, noutputs=12, seed=12345678, input = input)
  
  
  obs <- SIRdata
}


# Set up Libbi model ------------------------------------------------------

model<- libbi(SIRmodel, 
              nsamples = 1000, end_time = 73,
              nparticles = 1, obs = obs, 
              input = input, seed=1234,
              nthreads = 4,
              options = list(with="transform-initial-to-param"), verbose = verbose)


# Sample priors -----------------------------------------------------------

if (sample_priors) {
  prior <- sample(model, target = "prior")
}

# Run mcmc using the prior as the proposal --------------------------------

if (adapt_part || adapt_prop) {
  
  bi_prior <- sample(model, proposal="prior")
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
  adapted <- adapt_proposal(adapted, min=0.05, max=0.4, adapt = "both", scale = 1.1, max_iter = 3)
  
  get_block(adapted$model, "proposal_parameter")
}

if (save_results) {
  ModelTBBCGEngland::save_libbi(adapted, "rbi-test-adapted")
}

# Sample posterior using PMCMC --------------------------------------------

if (sample_post) {
  posterior <- sample(adapted, target = "posterior", nsamples = 1000, sample_obs = TRUE)
  
  plot(posterior)
}else{
  posterior <- adapted
}


# Sample posterior using SIR/SMC2 -----------------------------------------


if (use_sir_sampling) {
  posterior <- sample(posterior, target = "posterior", nsamples = 1000, sample_obs = TRUE,
                      options = list("sampler" = "sir", 
                                     "adapter" = "global",
                                     "tmoves" =  -1,
                                     "nmoves" = 1),
                      verbose = TRUE)
  
  p <- plot(posterior, plot = FALSE)
  
  p
}else{
  posterior <- adapted
}


# Predict states for all times. -------------------------------------------

if (pred_sample) {
  posterior <- predict(posterior, end_time = 120, noutputs = 120, verbose = FALSE)
  
  p <- plot(posterior, plot = FALSE)
  
  p
}

if (save_results) {
ModelTBBCGEngland::save_libbi(posterior, "rbi-test-model")
}
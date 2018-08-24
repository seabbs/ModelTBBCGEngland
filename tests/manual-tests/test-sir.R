library(ModelTBBCGEngland)


# Notes -------------------------------------------------------------------

## This function uses the SIR model supplied in the RBI package:
## https://github.com/sbfnk/RBi/blob/master/inst/SIR.bi
## All problem specific functionality has been disabled or ignored
## Results should be comparable to those found here: https://github.com/sbfnk/RBi.helpers/blob/master/vignettes/introduction.Rmd.rsp
## The primary aim of this script is to test the model generic aspects of the model fitting pipeline
# Fit model ---------------------------------------------------------------

model <- fit_model(model= "SIR", gen_data = TRUE, run_time = 16 * 7, plot_obs = FALSE,
                   sample_priors = FALSE, prior_samples = 1000, nparticles = 16, adaption_samples = 100, 
                   adapt_particles = TRUE, min_particles = 16, max_particles = 128,
                   adapt_proposal = TRUE, adapt_scale = 2, adapt_it = 10, adapt = "both",
                   fit = TRUE, posterior_samples = 500, thin = 1, burn_prop = 0, 
                   nthreads = parallel::detectCores(), verbose = TRUE, libbi_verbose = FALSE, fitting_verbose = FALSE,
                   save_output = FALSE, dir_path = "./tests/manual-tests/sir-results", dir_name = NULL)
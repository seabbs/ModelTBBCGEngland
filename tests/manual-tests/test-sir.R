library(ModelTBBCGEngland)


# Notes -------------------------------------------------------------------

## This function uses the SIR model supplied in the RBI package:
## https://github.com/sbfnk/RBi/blob/master/inst/SIR.bi
## All problem specific functionality has been disabled or ignored
## Results should be comparable to those found here: https://github.com/sbfnk/RBi.helpers/blob/master/vignettes/introduction.Rmd.rsp
## The primary aim of this script is to test the model generic aspects of the model fitting pipeline
# Fit model ---------------------------------------------------------------

model <- fit_model(model= "SIR", gen_data = TRUE, run_time = 16 * 7, plot_obs = FALSE,
                   sample_priors = TRUE, prior_samples = 1000, nparticles = 16,
                   adapt_particles = FALSE, min_particles = 16, max_particles = 256,
                   adapt_part_samples = 1000, adapt_part_it = 5, 
                   adapt_proposal = FALSE,  adapt_prop_samples = 1000, adapt_prop_it = 10,
                   adapt_scale = 2, adapt = "both",
                   fit = FALSE, posterior_samples = 10000, thin = 1, burn_prop = 0.25, 
                   nthreads = parallel::detectCores(), verbose = TRUE, libbi_verbose = FALSE, fitting_verbose = TRUE,
                   save_output = TRUE, dir_path = "./tests/manual-tests/sir-results", dir_name = NULL, seed = 12345678)


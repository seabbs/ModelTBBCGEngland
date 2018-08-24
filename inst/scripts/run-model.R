
# Install and load package ------------------------------------------------
library(ModelTBBCGEngland)


# Notes -------------------------------------------------------------------

# Fit model ---------------------------------------------------------------

model <- fit_model(model= "BaseLineModel", gen_data = TRUE, run_time = 73, plot_obs = TRUE,
                   sample_priors = FALSE, prior_samples = 1000, nparticles = 4, adaption_samples = 100, 
                   adapt_particles = TRUE, min_particles = 4, max_particles = 32,
                   adapt_proposal = TRUE, adapt_scale = 2, adapt_it = 10, adapt = "both",
                   fit = TRUE, posterior_samples = 500, thin = 1, burn_prop = 0, 
                   nthreads = parallel::detectCores(), verbose = TRUE, libbi_verbose = TRUE, fitting_verbose = TRUE, browse = TRUE,
                   save_output = TRUE, dir_path = "./vignettes/results/adapt", dir_name = NULL)


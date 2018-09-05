
# Install and load package ------------------------------------------------
library(ModelTBBCGEngland)


# Notes -------------------------------------------------------------------

# Fit model ---------------------------------------------------------------

model <- fit_model(model= "BaseLineModel", gen_data = FALSE, run_time = 73, plot_obs = TRUE,
                   sample_priors = TRUE, prior_samples = 1000, nparticles = 16,
                   adapt_particles = FALSE, min_particles = 4, max_particles = 128, adapt_part_samples = 1000, adapt_part_it = 3, 
                   adapt_proposal = TRUE,  adapt_prop_samples = 1000, adapt_prop_it = 5,
                   adapt_scale = 2, adapt = "shape",
                   fit = TRUE, posterior_samples = 10000, thin = 10, burn_prop = 0, 
                   nthreads = parallel::detectCores(), verbose = TRUE, libbi_verbose = FALSE, fitting_verbose = TRUE, browse = FALSE,
                   save_output = TRUE, dir_path = "./vignettes/results", dir_name = NULL)


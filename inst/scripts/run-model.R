
# Install and load package ------------------------------------------------
library(ModelTBBCGEngland)


# Notes -------------------------------------------------------------------

# Fit model ---------------------------------------------------------------

model <- fit_model(model= "BaseLineModel", gen_data = FALSE, run_time = 73, plot_obs = TRUE,
                   sample_priors = TRUE, prior_samples = 1000, nparticles = 4,
                   adapt_particles = TRUE, min_particles = 4, max_particles = 64, adapt_part_samples = 100, adapt_part_it = 5, 
                   adapt_proposal = TRUE,  adapt_prop_samples = 1000, adapt_prop_it = 10,
                   adapt_scale = 2, adapt = "both",
                   fit = TRUE, posterior_samples = 2000, thin = 1, burn_prop = 0, 
                   nthreads = parallel::detectCores(), verbose = TRUE, libbi_verbose = FALSE, fitting_verbose = TRUE, browse = FALSE,
                   save_output = TRUE, dir_path = "./vignettes/results/adapt", dir_name = NULL)


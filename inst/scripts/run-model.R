
# Install and load package ------------------------------------------------
library(ModelTBBCGEngland)


# Notes -------------------------------------------------------------------

# Fit model ---------------------------------------------------------------

model <- fit_model(model= "BaseLineModel", gen_data = FALSE, run_time = 73, plot_obs = TRUE,
                   sample_priors = TRUE, prior_samples = 1000, nparticles = 16,
                   adapt_particles = FALSE, min_particles = 4, max_particles = 128, adapt_part_samples = 1000, adapt_part_it = 3, 
                   adapt_proposal = TRUE,  adapt_prop_samples = 1000, adapt_prop_it = 3,
                   adapt_scale = 2, adapt = "both",
                   fit = TRUE, posterior_samples = 5000, thin = 1, burn_prop = 0, 
                   nthreads = parallel::detectCores(), verbose = TRUE, libbi_verbose = FALSE, fitting_verbose = TRUE, browse = FALSE,
                   save_output = TRUE, dir_path = "./vignettes/results/adapt", dir_name = NULL)


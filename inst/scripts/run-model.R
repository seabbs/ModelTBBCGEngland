
# Install and load package ------------------------------------------------
library(ModelTBBCGEngland)


# Notes -------------------------------------------------------------------

# Fit model ---------------------------------------------------------------

model <- fit_model(model= "BaseLineModel", gen_data = FALSE, run_time = 73, plot_obs = TRUE,
                   sample_priors = TRUE, prior_samples = 1000, nparticles = 1,
                   adapt_particles = FALSE, min_particles = 1, max_particles = 128, adapt_part_samples = 1000, adapt_part_it = 1, 
                   adapt_proposal = TRUE,  adapt_prop_samples = 1000, adapt_prop_it = 3,
                   adapt_scale = 4, adapt = "both",
                   fit = TRUE, posterior_samples = 25000, thin = 20, burn_prop = 0.2, 
                   scale_rate_treat = TRUE, 
                   years_of_age = NULL, noise = FALSE,
                   nthreads = parallel::detectCores(), verbose = TRUE, libbi_verbose = FALSE,
                   fitting_verbose = TRUE, browse = FALSE,
                   save_output = TRUE, dir_path = "./vignettes/results", dir_name = NULL)


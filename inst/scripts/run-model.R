library(ModelTBBCGEngland)


# Notes -------------------------------------------------------------------

# Fit model ---------------------------------------------------------------

model <- fit_model(model= "BaseLineModel", gen_data = FALSE, run_time = 73, plot_obs = TRUE,
                   sample_priors = TRUE, prior_samples = 1000, nparticles = 16, adaption_samples = 1000, 
                   adapt_particles = TRUE, min_particles = 16, max_particles = 256,
                   adapt_proposal = TRUE, adapt_scale = 1.1, adapt_it = 20, adapt = "size",
                   fit = FALSE, posterior_samples = 2000, thin = 1, burn_prop = 0, 
                   nthreads = parallel::detectCores(), verbose = TRUE, libbi_verbose = FALSE, fitting_verbose = FALSE,
                   save_output = TRUE, dir_path = "./vignettes/results/adapt", dir_name = NULL)


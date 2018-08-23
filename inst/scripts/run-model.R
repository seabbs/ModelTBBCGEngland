library(ModelTBBCGEngland)

# Fit model ---------------------------------------------------------------

model <- fit_model(model= "BaseLineModel", gen_data = FALSE, run_time = 73, time_scale = "year", plot_obs = TRUE,
                   sample_priors = TRUE, prior_samples = 1000, nparticles = 2, adaption_samples = 1000, adapt_particles = TRUE,
                   min_particles = 2, max_particles = 8, proposal_param_block = NULL, proposal_initial_block = NULL, 
                   adapt_proposal = TRUE, adapt_scale = 2, min_acc = 0.05, max_acc = 0.4, 
                   fit = FALSE, posterior_samples = 5000, thin = 0, burn_prop = 0, 
                   nthreads = parallel::detectCores(), verbose = TRUE, libbi_verbose = TRUE, browse = FALSE,
                   const_pop = FALSE, no_age = FALSE, no_disease = FALSE,
                   save_output = TRUE, dir_path = "./tests/manual-tests/results", dir_name = NULL)
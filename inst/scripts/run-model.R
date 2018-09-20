
# Install and load package ------------------------------------------------
library(ModelTBBCGEngland)


# Notes -------------------------------------------------------------------

# Fit model ---------------------------------------------------------------

model <- fit_model(model = "BaseLineModel", gen_data = FALSE, 
                   run_time = 73, time_scale = "year", plot_obs = TRUE,
                   sample_priors = TRUE, prior_samples = 10, nparticles = 4,
                   adapt_particles = FALSE, adapt_part_samples = 100, adapt_part_it = 1, 
                   min_particles = NULL, max_particles = NULL,
                   proposal_param_block = NULL, proposal_initial_block = NULL, 
                   adapt_proposal = FALSE, adapt_prop_samples = 100, adapt_prop_it = 1, adapt = "both",
                   adapt_scale = 2, min_acc = 0.05, max_acc = 0.4,
                   fit = FALSE, posterior_samples = 100, thin = 1, burn_prop = 0, sample_ess_at = 0.8,
                   rejuv_moves = 10, nthreads = 4, verbose = TRUE, libbi_verbose = TRUE, 
                   fitting_verbose = TRUE, pred_states = FALSE, browse = FALSE,
                   const_pop = FALSE, no_age = FALSE, no_disease = FALSE, scale_rate_treat = TRUE, years_of_age = NULL,
                   noise = TRUE,
                   save_output = FALSE, dir_path = "./vignettes/results", dir_name = NULL, reports = TRUE)


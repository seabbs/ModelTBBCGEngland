
# Install and load package ------------------------------------------------
library(ModelTBBCGEngland)


# Notes -------------------------------------------------------------------

# Fit model ---------------------------------------------------------------

model <- fit_model(model = "BaseLineModel", gen_data = FALSE, 
                   run_time = 73, time_scale = "year", plot_obs = TRUE,
                   sample_priors = TRUE, prior_samples = 1000, nparticles = NULL,
                   adapt_particles = FALSE, adapt_part_samples = 100, adapt_part_it = 1, 
                   min_particles = NULL, max_particles = NULL,
                   proposal_param_block = NULL, proposal_initial_block = NULL, 
                   adapt_proposal = FALSE, adapt_prop_samples = 250, adapt_prop_it = 4, adapt = "size",
                   adapt_scale = 1.2, min_acc = 0.05, max_acc = 0.3,
                   fit = FALSE, posterior_samples = 100, sample_ess_at = 0.2,
                   rejuv_moves = 5, nthreads = 4, verbose = TRUE, libbi_verbose = TRUE, 
                   fitting_verbose = TRUE, pred_states = FALSE, browse = FALSE,
                   const_pop = FALSE, no_age = FALSE, no_disease = FALSE, scale_rate_treat = TRUE, years_of_age = c(2000, 2004), age_groups = NULL, con_age_groups = c("childen"),
                   noise = TRUE,
                   save_output = FALSE, dir_path = "./vignettes/results", dir_name = NULL, reports = TRUE)


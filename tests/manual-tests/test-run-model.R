library(tictoc)
library(ModelTBBCGEngland)

tic()

tb_model <- fit_model(model= "BaseLineModel", gen_data = FALSE, run_time = 73, time_scale = "year", plot_obs = TRUE,
                      sample_priors = TRUE, prior_samples = 10, nparticles = 16, adaption_samples = 100,
                      adapt_particles = FALSE, adapt_proposal = FALSE, adapt_scale = 2,
                      min_acc = 0.05, max_acc = 0.4,  fit = TRUE, posterior_samples = 5000, nthreads = 4, 
                      verbose = TRUE, libbi_verbose = TRUE, browse = FALSE,
                      const_pop = FALSE, no_age = FALSE, no_disease = FALSE)


toc() 
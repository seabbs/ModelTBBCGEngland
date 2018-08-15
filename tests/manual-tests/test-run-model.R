library(tictoc)
library(ModelTBBCGEngland)

tic()

tb_model <- test_model(model= "BaseLineModel", gen_data = FALSE, run_time = 74, time_scale = "year", plot_input_data = FALSE,
                      sample_priors = TRUE, prior_samples = 1000, nparticles = NULL, adapt_particles = FALSE,
                      adapt_proposal = FALSE, min_acc = 0.05, max_acc = 0.4, fit = FALSE, posterior_samples = 100,
                      nthreads = 4, verbose = TRUE, libbi_verbose = FALSE, browse = FALSE,
                      const_pop = FALSE, no_age = FALSE, no_disease = FALSE)


toc()
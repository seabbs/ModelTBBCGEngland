library(tictoc)
library(ModelTBBCGEngland)

tic()

tb_model <- test_model(model= "BaseLineModel", gen_data = TRUE, run_time = 73, time_scale = "year", plot_input_data = TRUE,
                      sample_priors = TRUE, prior_samples = 1000, nparticles = 16, adapt_particles = TRUE,
                      adapt_proposal = TRUE, min_acc = 0.05, max_acc = 0.4, fit = TRUE, posterior_samples = 100,
                      nthreads = 4, verbose = TRUE, libbi_verbose = FALSE, browse = FALSE,
                      const_pop = FALSE, no_age = FALSE, no_disease = FALSE)


toc()
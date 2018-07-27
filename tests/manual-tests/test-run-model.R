library(tictoc)
library(ModelTBBCGEngland)

tic()

tb_model <- test_model(model= "BaseLineModel", gen_data = TRUE, run_time = 10, time_scale = "year", plot_input_data = TRUE,
                      sample_priors = TRUE, prior_samples = 10, nparticles = NULL, adapt_particles = FALSE,
                      adapt_proposal = FALSE, min_acc = 0.05, max_acc = 0.4, fit = FALSE, posterior_samples = 100,
                      nthreads = 4, verbose = TRUE, libbi_verbose = TRUE, browse = FALSE)


toc()
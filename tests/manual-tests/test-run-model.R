library(tictoc)
library(ModelTBBCGEngland)

tic()

tb_model <- test_model(model= "BaseLineModel", gen_data = TRUE, run_time = 74, plot_input_data = TRUE,
                      sample_priors = FALSE, nsamples = 1000, nparticles = 1000, adapt_particles = FALSE,
                      adapt_proposal = FALSE, min_acc = 0.05, max_acc = 0.4, fit = FALSE,
                      nthreads = 4, verbose = TRUE, libbi_verbose = FALSE)


toc()
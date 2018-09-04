library(ModelTBBCGEngland)
library(tidyverse)

prior_params <- readRDS("vignettes/results/adapt/model-run-2018-09-04_11:02:42/data/prior-params.rds")
libbi_model <- ModelTBBCGEngland::read_libbi("vignettes/results/adapt/model-run-2018-09-04_11:02:42/libbi/posterior.rds")

plot_param(libbi_model, prior_params = prior_params, scales = "free")


plot_state(libbi_model, c("YearlyInc", "YearlyHistPInc"), burn_in = 50)

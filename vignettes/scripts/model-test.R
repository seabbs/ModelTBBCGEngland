
# Load packages -----------------------------------------------------------

library(tidyverse)
library(rbi)
library(rbi.helpers)


# Load the model ----------------------------------------------------------

model_file <- file.path(getwd(), "inst/bi/BaseLineModel.bi")

tb_model <- bi_model(model_file)


# Generate data from the model --------------------------------------------

tb_data <- bi_generate_dataset(tb_model, end_time = 12 * 75, noutputs = 16, seed = 12345678)

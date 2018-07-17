
# Load packages -----------------------------------------------------------

library(tidyverse)
library(rbi)
library(rbi.helpers)


# Load the model ----------------------------------------------------------

model_file <- system.file(package="ModelTBBCGEngland", "BaseLineModel.bi")

tb_model <- bi_model(model_file)

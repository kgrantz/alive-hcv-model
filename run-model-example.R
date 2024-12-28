
# code to run example transmission model
# - Longitudinal trajectories
# - Age-specific mixing
# - post-tx behavioral changes
# - no class-specific treatment

library(tidyverse)
library(lcmm)

source("indiv-model-functions.R")

load("example_params.Rdata")

init_model(pop_dir = "populations",
           out_dir = "output/",
           betas = c(t0=0.1, t5474=0.05, t9128=0.03),
           tx_rates = c(t0=0, t9128=0.1, t10958=0.2),
           mixing_df = mij,
           start_time = 0,
           end_time = nt,
           posttx_effects = TRUE,
           long_traj = TRUE,
           class_tx = NULL)

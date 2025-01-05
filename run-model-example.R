
# code to run example transmission model
# - Longitudinal trajectories
# - Age-specific mixing
# - post-tx behavioral changes
# - no class-specific treatment

library(tidyverse)
library(lcmm)

source("indiv-model-functions.R")

load("example_params.Rdata")

# this will run a model from scratch (i.e., starting from the original population
# object and time = 0).
total_pop <- read_parquet("populations/total_pop.parquet")
load("populations/total_pop_ids.Rdata")

init_model(total_pop = total_pop,
           out_id = out_id,
           out_dir = "output/scn1_initial",
           betas = c(t0=7.5e-5, t9135=6e-5), 
           tx_rates = c(t0=0), # no treatment prior to 2014
           mixing_df = mij,
           start_time = 0,
           end_time = 9135, # for running through end of 2014
           posttx_effects = TRUE,
           long_traj = TRUE,
           class_tx = NULL,
           beta_fit = FALSE,
           summary_out = FALSE)

# this code will run model beginning from the end of the prior simulation
# this is useful, e.g., for initiating different treatment strategies at the
# same time point for equivalent burn-in/calibration periods
init_model(pop_dir = "output/scn1_initial",
           out_dir = "output/scn1_tx1",
           betas = c(t0=7.5e-5, t9135=6e-5), 
           tx_rates = c(t0=0, t9135=0.1, t10962=0.3),
           mixing_df = mij,
           start_time = 9135,
           end_time = nt,
           posttx_effects = TRUE,
           long_traj = TRUE,
           class_tx = NULL,
           beta_fit = FALSE,
           summary_out = TRUE)

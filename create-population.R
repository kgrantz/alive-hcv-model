
## code to create population sets for use in HCV modeling

library(tidyverse)
library(lcmm)

source("indiv-model-functions.R")

load("example_params.Rdata")
mu <- mu1

out_dir <- "populations/"
stopifnot(dir.exists(out_dir))

# just for testing, make a small number of small populations
# this is overwriting values in params
nsim = 5 
n0 = 100

# format LCMM object for making individual parameters
mod_risk <- format_lcmm(mspl_ci)

# loop across number of simulations (population sets) to run
for(i in 1:nsim){
  
  # make population
  pop <- make_population(n0=n0) %>%
    mutate(nsim = i)
  
  # add individual parameters
  pop <- make_indiv_params(pop, 
                           alpha1=sc_alpha, 
                           beta1=sc_beta, 
                           alpha2=txred_alpha, 
                           beta2=txred_beta, 
                           risk_class=mod_risk)
  
  # save output
  save(pop, file=glue("{out_dir}/pop_{i}.Rdata"))
  
}
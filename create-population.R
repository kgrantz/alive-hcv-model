
## code to create population sets for use in HCV modeling
##
## starts by creating one large set of individuals,
## and then creating unique lists of subsets of IDs to read in
## for individual simulations

library(tidyverse)
library(lcmm)

source("indiv-model-functions.R")

load("example_params.Rdata")

# must have a directory to store output - recommend using `populations`
out_dir <- "populations/"
stopifnot(dir.exists(out_dir))

# format LCMM object for making individual parameters
mod_risk <- format_lcmm(mspl_rf)

# format death rates
mu_format = format_mu(mu)

## -----------------------------------------------------------
## step 1: make set of individuals to draw for population sets
## -----------------------------------------------------------

mult_fac = 5 # number of populations of size n0 to combine to make set of individuals; recommend this be >= nsim/10
total_pop_size = n0*mult_fac 
pop_iter = mult_fac # build the total pop in several iterations, so smooths out some of the variability in total population size
iter_pop_size = total_pop_size/pop_iter
total_pop = data.frame()
final_pop_size = rep(0, pop_iter)

for(i in 1:pop_iter){
  
  print(glue("iteration {i} startion {date()}"))
  
  # make population
  tmp = make_population(n0 = iter_pop_size)
  
  # adjust ID for number of individuals already in total pop
  if(nrow(total_pop)>0){
    tmp$id = tmp$id + max(total_pop$id)
  }
  
  final_pop_size[i] = length(unique(tmp$id))
  
  # combine into new pop
  total_pop = bind_rows(total_pop, tmp)
}

# add individual parameters
total_pop <- make_indiv_params(total_pop, 
                               alpha1=sc_alpha, 
                               beta1=sc_beta, 
                               alpha2=txred_alpha, 
                               beta2=txred_beta, 
                               sct1=sctime_min,
                               sct2=sctime_max,
                               risk_class=mod_risk)

# save total_pop file
write_parquet(total_pop, glue("{out_dir}/total_pop.parquet"))
save(final_pop_size, file=glue("{out_dir}/final_pop_size.Rdata"))

## -----------------------------------------------------------
## step 2: draw IDs for each population set
## -----------------------------------------------------------


total_pop_ids <- unique(total_pop$id)
pop_size = sample(final_pop_size, nsim, replace=TRUE)

out_id <- list()

for(i in 1:nsim){
  out_id[[i]] = sample(total_pop_ids, pop_size[i], replace=T)  
}

save(out_id, file=glue("{out_dir}/total_pop_ids.Rdata"))


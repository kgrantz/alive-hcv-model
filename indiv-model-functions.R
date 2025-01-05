
library(tidyr)
library(dplyr)
library(tibble)
library(readr)
library(lubridate)
library(glue)
library(lcmm)
library(truncnorm)
library(arrow)

## ................................................
## DEMOGRAPHY FUNCTIONS........................####
## ................................................


#' Function to create individuals in population generation
#' 
#' @param n numeric, number of individuals to create
#' @param t numeric, time at which creation is occuring
#' @param starting_id numeric, the first ID that should be assigned to a newly created person. Should be max(id)+1 of existing population
#' @param max_time numeric, maximum time point
#' @param time_step numeric, length of time steps (default 7d)
#' @param class_probs list, with $levels (1, 2, 3) and $probs (proportion of pop to be in each class)
#' @param age_class list of length 3 with empirical distributions of age of enrolment for each class
#' @param age_idu_class list of length 3 with empirical distributions of age at first IDU for each class
#' @param hcv_class vector of length 3 with SCALED probability of HCV seropositivity at enrolment for each class; scaled given the curve dictated by hcv_age_prm
#' @param hcv_vir numeric, probability of being viremic given seropositive at enrolment for entire cohort
#' @param hcv_age_prm numeric, parameter controlling age-specific prevalence curve at baseline
#' @param age_brks age breaks to use in creating age categories
#' @param death_df formatted df of deaths rates, output of format_mu
#' 
make_individuals <- function(n, t, starting_id, 
                             max_time=nt, 
                             time_step=dt, 
                             class_probs=classes, 
                             age_class=age_enrol, 
                             age_idu_class=age_idu, 
                             hcv_class=hcv_prev_scaled, 
                             hcv_vir_pr=hcv_vir_prob,
                             hcv_age_prm=hcv_age_param,
                             age_brks=age_breaks, 
                             death_df = mu_format){
  
  if(n>0){
    # create df of just IDs to start
    xx <- data.frame(
      id = starting_id:(starting_id+n-1)  
    )
    
    # assign class membership, age, age first IDU, baseline HCV prev
    nsamp <- nrow(xx)
    xx$class <- sample(class_probs$levels, nrow(xx), prob=class_probs$probs, replace=T)
    xx$age_enrol[xx$class==1] <- sample(age_class[[1]], sum(xx$class==1), replace=T)
    xx$age_enrol[xx$class==2] <- sample(age_class[[2]], sum(xx$class==2), replace=T)
    xx$age_enrol[xx$class==3] <- sample(age_class[[3]], sum(xx$class==3), replace=T)
    xx$age_firstidu[xx$class==1] <- sample(age_idu_class[[1]], sum(xx$class==1), replace=T)
    xx$age_firstidu[xx$class==2] <- sample(age_idu_class[[2]], sum(xx$class==2), replace=T)
    xx$age_firstidu[xx$class==3] <- sample(age_idu_class[[3]], sum(xx$class==3), replace=T)
    age_class_num = as.numeric(cut(xx$age_enrol, age_breaks, right=FALSE))
    
    # assign baseline HCV prev + viremia
    # if age_firstidu <= age_enrol, draw HCV prev as 'migration'
    # if age_firstidu > age_enrol, make HCV negative 
    # hcv prevalence is adjusted by an age specific curve (estimated using empirical age specific prevalence)
    xx$hcv_ab <- 0L
    ind <- which(xx$age_firstidu <= xx$age_enrol)
    hcv_ab_pr = (1-exp(-hcv_age_prm*age_class_num)) * hcv_class[xx$class]
    hcv_ab_pr[hcv_ab_pr>0.99] = 0.99
    xx$hcv_ab[ind] = rbinom(length(ind), 1, hcv_ab_pr[ind])
    xx$hcv_vir = 0L
    ind <- which(xx$hcv_ab==1)
    xx$hcv_vir[ind] = rbinom(length(ind), 1, hcv_vir_pr)
    
    # expand out IDs + times and merge in the
    # individual ID characteristics
    tmp <- 
      expand.grid(
        id = xx$id,
        time = seq(t, max_time, time_step)
      ) 
    tmp <- left_join(tmp, xx, by="id")
    
    # calculate age/age_cat over time
    tmp$age <- tmp$age_enrol + (tmp$time-t)/365.25
    tmp$age_cat <- cut(tmp$age, age_brks, right=FALSE)
    
    # filter so that individuals are only included when age_firstidu<=age
    # that is, entry when age is first >=age_firstidu
    tmp <- tmp[tmp$age>=tmp$age_firstidu,]
    
    # add death rate + probability
    tmp <- left_join(tmp, death_df, by=c("time", "age_cat"))
    
    # for each individual, draw death events
    # filter once reach first death
    tmp$death <- rbinom(nrow(tmp), 1, tmp$pr_death)
    tmp <- tmp %>%
      group_by(id) %>%
      filter(cumsum(death)<1) %>%
      # create tau, time since entry
      mutate(tau = time - min(time))
    
  }else{
    tmp <- NULL
  }
  return(tmp)
}


#' Function to create population, using make_individuals
#' assumes environment has default parameters of make_individuals loaded
#' (from params/params-scenario2.Rdata)
#' @param n0 number of people to create at time 0
#' 

make_population <- function(n0, max_time=nt, time_step=dt){
  
  # start with baseline pop at t=0
  pop <- make_individuals(n=n0, t=0, starting_id=1)
  
  # make the times vector
  times <- seq(0, max_time, time_step)[-1]
  
  # for each time step, 
  # draw number of births and pass
  # arguments to make_individuals
  for(i in times){
    tmp <- pop[pop$time==i,]
    starting_id = max(pop$id)+1
    n = sum(rbinom(nrow(tmp), 1, tmp$pr_death))
    t = i
    add <- make_individuals(n=n, t=t, starting_id=starting_id)
    pop <- bind_rows(pop, add)
  }
  
  return(pop)
}

#' convenience function to reformat LCMM predictions for
#' merging with population
#' @param mod LCMM model output
#' 

format_lcmm <- function(mod){
  
  out <- mod$pred %>%
    as_tibble() %>%
    # add time_since_enrol = tau in model notation
    bind_cols(mod$times) %>%
    rename(tau = time_since_enrol) %>%
    pivot_longer(-tau, names_sep="_", names_to=c("var", "quant", "class")) %>%
    mutate(quant = paste0("quant", quant),
           class = as.character(parse_number(class))) %>%
    pivot_wider(names_from = quant) %>%
    # make empirical std deviation
    mutate(sd_gamma = (quant97.5 - quant2.5)/1.96) %>%
    rename(mean_gamma=quant50) %>%
    # use only first 13 years
    filter(tau <= 13*365.25) %>%
    # retain only tau, class, mean, sd to feed into rnorm later
    select(tau, class, mean_gamma, sd_gamma)
  
  return(out)  
}

#' convenience function to reformat mortality rates for
#' merging with population
#' @param mu df, age categories and age-specific death rates
#' @param max_time numeric, maximum time point
#' @param time_step numeric, length of time steps (default 7d)
#' @param year_t0 calendar year for t=0 (default 1990)
#' @param year_brks year breaks to make year categories for merging death rates
#' @param year_lbl year category labels for convenience

format_mu <- function(mu,
                      max_time=nt, 
                      time_step=dt, 
                      year_t0=year0, 
                      year_brks=year_breaks, 
                      year_lbl=year_labels){
  
  # determine if mu contains variable year_cat
  if(is.null(mu$year_cat)){
    # if it doesn't, just expand out age cats/time steps
    out <- expand.grid(
      time = seq(0, max_time, time_step), 
      age_cat = unique(mu$age_cat)
    ) %>%   
    left_join(mu %>% select(age_cat, mu), by=c("age_cat"))    
  }else{
    # if year categories do exists, create year caategories for time steps
    # and merge in mu this way
    out <- expand.grid(
        time = seq(0, max_time, time_step), 
        age_cat = unique(mu$age_cat)
      ) %>% 
      # make year/year category for each time step
      mutate(year = floor(year0+time/365.25), 
             year = if_else(year>2018, 2018, year),  
             year_cat = cut(year, breaks=year_brks, labels=year_lbl, right=TRUE)) %>%
      left_join(mu %>% select(year_cat, age_cat, mu), by=c("age_cat", "year_cat")) 
  }

  out <- out %>%
    # calculate pr_death
    mutate(pr_death = 1 - exp(-time_step*mu/365.25)) %>%
    # keep only time, age_cat, pr_death for merging into pop
    select(time, age_cat, pr_death)   
  
  return(out)  
}

#' function to add individual parameters to population
#' including spontaneous clearance probability, post-tx risk reductions,
#' time to progress to chronic infection, and risk trajectories
#' @param pop df, population created by make_population
#' @param alpha1 shape param 1 of beta distribution of spontaneous clearance probability 
#' @param beta1 shape param 2 of beta distribution of spontaneous clearance probability
#' @param alpha2 shape param 1 of beta distribution of post-tx risk reduction
#' @param beta2 shape param 2 of beta distribution of post-tx risk reduction
#' @param risk_class output of format_lcmm, for making risk trajectories
#' 
make_indiv_params <- function(pop, 
                              delta_t=dt,
                              alpha1 = sc_alpha, beta1 = sc_beta,
                              alpha2 = txred_alpha, beta2 = txred_beta,
                              sct1 = sctime_min, sct2=sctime_max,
                              risk_class = mod_risk){
  
  xx <- data.frame(id=unique(pop$id))
  nsamp = nrow(xx)
  xx$pr_clearance_total <- rbeta(nsamp, alpha1, beta1)
  xx$posttx_red <- rbeta(nsamp, alpha2, beta2)
  xx$t_become_chronic <- sample(sct1:sct2, nsamp, replace=T)
  
  # recall pr_clearance_total = per-person rate over t_become_chronic
  # therefore rate = pr_clearance_total/t_become_chronic
  # and the probability = 1 - exp(-dt*pr_clearance_total/t_become)
  xx$pr_clearance = 1 - exp(-delta_t*xx$pr_clearance_total/xx$t_become_chronic)
  
  pop <- left_join(pop, xx, by="id")
  
  # add risk class trajectories (per ID and time point; make new tau so
  # that any tau > 13y = 13y, since trajectories only extend to 13 years
  # post IDU initiation)
  pop$tau_orig = pop$tau
  pop$tau[pop$tau>13*365.25] = floor(13*365.25)
  
  xx <- left_join(pop[,c("id","tau","class")], risk_class, by=c("tau", "class"))
  pop$gamma <- rtruncnorm(nrow(pop), mean=xx$mean_gamma, sd=xx$sd_gamma, a=0)

  return(pop)
}

#' function to make a population set from total set of individuals
#' and vector of IDs to be included in the population
#' re-orders IDs to give a unique identifier for each population
#' 
#' @param total_pop data frame of total set of individuals from which to boostrap populations
#' @param ids vector of sampled IDs from total_pop to use to create single population iteration
make_pop_from_ids <-  function(total_pop, ids){
  tmp = data.frame(id = ids) %>%
    mutate(id_new = 1:n()) %>%
    left_join(total_pop, by="id") %>%
    fill(id_new) %>%
    select(-id) %>%
    rename(id = id_new)
  return(tmp)
}

## ................................................
## TRANSMISSION FUNCTIONS......................####
## ................................................


#' function to compute and add treatment rates to pop df
#' @param pop df, population created by make_population and make_indiv_params
#' @param tx_rates named vector of annual per-person treatment rates; assumes step function
#' betweeen time points. e.g., tx_rates=c(t0=0.2, t365=0.05) means annual per-person treatment rate 
#' will be 0.2 from t0 to t364, and then 0.05 from t365 onwards.
#' tx_rates must have at least one entry for t0, all other times flexible.
#' In applications for ALIVE, treatment should be 0 at t0 (1990), and then a low rate
#' (~ 0.1) from 2015-2019 (start 2015: t=9142) and then set to whatever desired
#' levels to model after 2020 (t=10962)
#' 
make_tx_rates <- function(pop, tx_rates){
  
  # check that there is a t0 in tx_rates
  if(is.na(tx_rates["t0"])){
    stop("`tx_rates` must contain an entry for t0")
  }
  
  # create vector of tx_rates:
  tx <- tibble(
    time = parse_number(names(tx_rates)),
    tx_rate = tx_rates
  ) %>%
    # add in all times from population
    bind_rows(tibble(time = unique(pop$time))) %>%
    # remove duplicates of time 
    distinct(time, .keep_all=TRUE) %>%
    # rearrange in time order
    arrange(time) %>%
    # fill treatment rates down
    fill(tx_rate)
  
  # add the treatment rate into pop df 
  pop <- left_join(pop, tx, by="time")
  
  return(pop)
}

#' function to compute time-specific beta prefactors using linear
#' extrapolation between provided timepoints
#' @param pop df, population created by make_population and make_indiv_params
#' @param betas named vector of beta prefactors; names should be the starting
#' time point of linear extrapolations. In practical application in fitting, the
#' time points to use are t0 (1990), t5488 (start 2005) and t9135 (last time point of 
#' 2014, before 2015 tx starts)
make_beta <- function(pop, betas){
  
  # check that there is a t0 in betas
  if(is.na(betas["t0"])){
    stop("`betas` must contain an entry for t0")
  }
  
  # pull out the times where beta changes
  beta_times = parse_number(names(betas))
  
  # create vector of tx_rates:
  bt <- tibble(
    time = beta_times,
    beta = betas
  )  
  
  # if only one beta value is given, make beta fixed for all time points
  if(nrow(bt)==1){
    pop$beta <- bt$beta
  }else{
    # otherwise, calculate slopes between time points
    # beta will remain constant (i.e., slope=0) after last time point in beta
    bt <- bt %>%
      mutate(slope = c(diff(beta)/diff(time), 0)) %>%
      # add in all times from population
      bind_rows(tibble(time = unique(pop$time))) %>%
      # remove duplicates of time 
      distinct(time, .keep_all=TRUE) %>%
      # rearrange in time order
      arrange(time) %>%
      # fill beta down
      fill(beta, slope) %>%
      # create categories for the times to get x
      mutate(cat = cut(time, beta_times, right=F)) %>%
      group_by(cat) %>%
      mutate(x = time - min(time),
             beta = beta + slope*x) %>%
      ungroup %>%
      # keep just time and slope
      select(time, beta)
    
    # merge with pop
    pop <- pop %>%
      left_join(bt, by="time")
  }
  
  return(pop)
}


#' function to calculate treatment probabilities and draw 
#' treatment events per person and per time step in model transmission
#' @param tmp df, population object at one time step in model transmission
#' @param class_tx_prob vector of length 3 with class-specific relative rates
#' of treatment; if NULL, all classes have same treatment rate
#' 
make_treatment <- function(tmp, class_tx_prob=NULL){
  
  # determine if any treatment is occurring
  # if so, run treatment module
  if(sum(tmp[,"tx_rate"],na.rm=T)>0){
    
    # first, define the class-specific treatment rates, 
    # if class_tx_prob exists
    if(!is.null(class_tx_prob)){
      
      # balance the treatment rates
      # first find proportion in each class
      prop_class <- table(tmp[,"class"])/nrow(tmp) # tmp %>% count(class) %>% mutate(p=n/sum(n))

      # calculate the treatment scalar
      tx_scalar = 1 / sum(class_tx_prob * prop_class) #prop_class$p)
      
      # make tx_rate_scaled
      tx_rate_scaled = tmp[,"tx_rate"] * tx_scalar * class_tx_prob[tmp[,"class"]]
      
    }else{
      
      # if no class scaling, just make tx_rate_scaled = tx_rate
      tx_rate_scaled = tmp[,"tx_rate"]
      
    }
    
    # create probability of treatment + draw treatment instance + treatment success
    pr_tx = 1 - exp(-dt*tx_rate_scaled/365.25)
    tmp_tx = rep(0L, nrow(tmp))
    ind = which(tmp[,"st_char_prev"]==3)
    tmp_tx[ind] = rbinom(length(ind), 1, pr_tx[ind])
    tmp_tx_success = rep(0L, nrow(tmp))
    ind = which(tmp_tx==1)
    tmp_tx_success[ind] = rbinom(length(ind), 1, tx_eff)
    
  }else{
    # if treatment rates are all 0, then pr_tx=tmp_tx = 0, tmp_tx_success=NA for all
    pr_tx = 0
    tmp_tx = 0
    tmp_tx_success = 0
  }
  
  tmp <- cbind(tmp, pr_tx, tmp_tx, tmp_tx_success)
  
  return(tmp)
}


#' function to calculate infection probabilities and draw
#' infection events per person and per time step in model transmission
#' @param tmp df, population object at one time step in model transmission
#' @param mixing_df df of age class i, j, and mixing intensity between ages
#' 
#' How to calculate FOI:
#' - all acute + chronic individuals contribute to FOI
#' - all individuals can be infected
#' - FOI for individual i in age cat a_i:
#'    \lambda = \beta * \gamma_i * \sum_{j in Infected} \gamma_j * \K_{ai,aj}
#' - if we break down infected into Infected_aj=1 (i.e., single age categories), you get:
#'    \lambda = \beta * \gamma_i * \sum_{aj cats} \sum_{j in Infected_aj} \gamma_j * \K_{ai,aj=aj}
#' - for a given individual of age ai, K_{ai,aj=aj} is a constant:
#'    \lambda = \beta * \gamma_i * \sum_{aj cats} \K_{ai,aj=aj} \sum_{j in Infected_aj} \gamma_j
#' - therefore, to get FOI term for a given individual, we can:
#'   - sum gamma for all infectious (A or C) individuals in each category (inj_age_j)
#'   - for each age_cat_i, calculate sum(K_ai=i,aj * inf_age_j), i.e. summed contribution of
#'     all j cats to each i cat's FOI (tmp_foi)
#'   - merge these terms into tmp by age_cat_i
#'   - if the individual is infected, subtract gamma_i*K_ai,ai (so a person doesn't contribute to
#'   their own FOI)
#'   - multiple FOI term * beta * gamma_i * (1-epsilon*posttx_red)
#' - Draw infection events based on lamba rate
#' 
make_infection <- function(tmp, mixing_df, assort_mix){
  
  # If there is no class-assortative mixing, just calculate by age class:
  if(is.null(assort_mix)){
      
    # calculate sum(gamma) for each aj
    inf_ind <- which(tmp[,"st_char_prev"] > 1)
    inf_agej <- tapply(tmp[inf_ind, "gamma"]*(1-tmp[inf_ind,"oc_epsilon"]*tmp[inf_ind,"posttx_red"]), tmp[inf_ind, "age_cat"], sum)
    
    # calculate sum of Kai,aj*sum_aj(gamma) for each ai
    sum_gamma = inf_agej[match(as.numeric(mixing_df$agep_cat), names(inf_agej))] 
    foi_sum_age = tapply(sum_gamma*mixing_df$mij, as.numeric(mixing_df$agec_cat), function(x)(sum(x,na.rm=T)))
  
    # merge this sum into tmp by age_cat_i
    foi_sum <- foi_sum_age[match(tmp[,"age_cat"], names(foi_sum_age))]
  
  }
  
  # If there is class-assortative mixing, calculate by age and risk class
  if(!is.null(assort_mix)){
    
    # make age cat/class index
    tmp_age_class = interaction(tmp[,"age_cat"], tmp[,"class"])
    
    # calculate sum(gamma) for each aj, class
    inf_ind <- which(tmp[,"st_char_prev"] > 1)
    
    inf_agej <- tapply(
      tmp[inf_ind, "gamma"]*(1-tmp[inf_ind,"oc_epsilon"]*tmp[inf_ind,"posttx_red"]), 
      tmp_age_class[inf_ind],
      sum
    )
    
    # calculate sum of Kai,aj*sum_aj(gamma) for each ai, s
    sum_gamma = inf_agej[match(mixing_df$agep_cat_class, names(inf_agej))] 
    foi_sum_age = tapply(sum_gamma*mixing_df$mij, mixing_df$agec_cat_class, function(x)(sum(x,na.rm=T)))
    
    # merge this sum into tmp by age_cat_i, class
    foi_sum <- foi_sum_age[match(tmp_age_class, names(foi_sum_age))]
  }
  
  # calculate FOI for each person
  # if an individual is infected, remove their contribution to FOI
  # otherwise foi_sum remains unchanged
  foi_sum[inf_ind] = foi_sum[inf_ind] - tmp[inf_ind, "indiv_mij"]*tmp[inf_ind, "gamma"]*(1-tmp[inf_ind,"oc_epsilon"]*tmp[inf_ind,"posttx_red"])
  lambda = tmp[,"beta"]*tmp[,"gamma"]*foi_sum*(1-tmp[,"oc_epsilon"]*tmp[,"posttx_red"])
  pr_infection = 1 - exp(-dt * lambda)
  tmp_infection = rbinom(nrow(tmp), 1, pr_infection)
  
  tmp <- cbind(tmp, tmp_infection)
  return(tmp)
}


#' function to update individual infections states based on model events
#' @param df df, population object at one time step in model transmission
#' @param st_char_prev column in df containing the previous time point's state (N, A, C)
#' @param tmp_tx column in df, binary 0/1 whether individual was treated
#' @param tmp_tx_success column in df, binary 0/1 whether treatment was successful
#' @param tmp_chronic column in df, binary 0/1 whether individual progressed from A to C
#' @param tmp_clearance column in df, binary 0/1 whether individual cleared infection
#' @param tmp_infection column in df, binary 0/1 whether individual acquired infection
#' @param oc_epsilon, column in df, binary 0/1 whether individual has ever been treated
#' 
#' How to update states:
#' **starting state N:**
#'  - no infection: N
#'  - infection: A (t_since_infection=0)
#'**starting state A:**
#'  - clearance + infection: A (t_since_infection=0)
#'  - clearance + no infection: N
#'  - no clearance + chronic + inf/no inf: C
#'  - no clearance + no chronic + inf/no inf: A (t_since_infection +1)
#' **starting state C:**
#'  - treatment + successful + inf: A (t_since_infection 0)
#'  - treatment + successful + no inf: N
#'  - treatment + failure + inf/no inf: C
#'  - no treatment + inf/no inf: C
#' note that this schema functionally ranks clearance ahead of
#' progression to chronic infection (i.e., don't move to C if
#' you also clear in the same time step)
#' also says that any infection occurs after tx/clearance as reinfection,
#' not prior to
#' 
update_state <- function(df){
  
  # update infection outcomes:
  df[df[,"tmp_infection"]==1, "oc_sero_pos"] = 1 # if infected, oc_sero_pos should become 1; otherwise keep previous value
  df[,"oc_n_inf"] = df[,"oc_n_inf"] + df[,"tmp_infection"]

  ind <- which(df[,"st_char_prev"]==1)
  df[ind, "oc_n_new_inf"] = df[ind,"oc_n_new_inf"] + df[ind,"tmp_infection"]
  
  ind <- which(df[,"st_char_prev"]>1)
  df[ind,"oc_n_superinf"] = df[ind,"oc_n_superinf"] + df[ind,"tmp_infection"]
  
  df[,"oc_n_clearance"] = df[,"oc_n_clearance"] + df[,"tmp_clearance"]
  
  # update treatment outcomes
  df[,"oc_n_tx"] = df[,"oc_n_tx"] + df[,"tmp_tx"]
  df[,"oc_n_tx_success"] = df[,"oc_n_tx_success"] + df[,"tmp_tx"]*df[,"tmp_tx_success"]
  df[,"oc_epsilon"] = as.numeric(df[,"oc_n_tx"]>0)
  df[,"oc_n_posttx_inf"] = df[,"oc_n_posttx_inf"]+df[,"tmp_infection"]*df[,"oc_epsilon"]
  ind <- which(df[,"st_char_prev"]==1)
  df[ind,"oc_n_posttx_new_inf"] = df[ind,"oc_n_posttx_new_inf"]+df[ind,"tmp_infection"]*df[ind,"oc_epsilon"]
  
  # update infection states
  df[,"st_char"] = df[,"st_char_prev"] # encompasses N>N (no infection), A>A (no clearance+no chronic+inf/no inf, clearance+inf+chronic/no chronic), C>C (no treatment, treatment failure)
  df[df[,"st_char_prev"]==1 & df[,"tmp_infection"]==1,"st_char"] = 2 # new infection
  df[df[,"st_char_prev"]==2 & df[,"tmp_clearance"]==1 & df[,"tmp_infection"]==0,"st_char"] = 1 # clearing + not infected
  df[df[,"st_char_prev"]==2 & df[,"tmp_clearance"]==0 & df[,"tmp_chronic"]==1,"st_char"] = 3 # no clearance + becoming chronic
  df[df[,"st_char_prev"]==3 & df[,"tmp_tx"]==1 & df[,"tmp_tx_success"]==1 & df[,"tmp_infection"]==1,"st_char"] = 2 # treated + successful + reinfected
  df[df[,"st_char_prev"]==3 & df[,"tmp_tx"]==1 & df[,"tmp_tx_success"]==1 & df[,"tmp_infection"]==0,"st_char"] = 1 # treated successfully + no infection

  # update time since infection for acute individuals
  # time_since infection should be 0 for basically everyone,
  # except those who stayed in A and did not clear/get reinfected
  ind = which(df[,"st_char_prev"]==2 & df[,"st_char"]==2 & df[,"tmp_clearance"]==0)
  df[ind,"t_since_infection"] = df[ind,"t_since_infection"] + dt
  df[-ind, "t_since_infection"] = 0
  
  return(df)

}


#' function to initialize a model run across all populations
#' @param pop_dir directory where individual population RData files are stored - to be used if resuming from other model output
#' @param total_pop df of total population, to be used if building population sets within init_model
#' @param out_id list of IDs to be drawn from total_pop, to be used if building population sets within init_model
#' @param out_dir directory where simulation output will go
#' @param betas named vector of beta prefactors, to be passed to make_beta function
#' @param tx_rates named vector of treatment rates, to be passed to make_tx_rates function
#' @param mixing_df df of age class i, j, and mixing intensity between ages, to be passed to make_infection function
#' @param start_time time point to start running the model; if >0, assumed that pop_dir points to prior model output.
#' start_time should be the same as end_time for the output from which you are resuming the run (i.e., t=9135)
#' @param end_time time point to stop the model run, useful for calibrating betas (stop at t=9135)
#' @param posttx_effects logical for whether to include post-tx reductions in transmission
#' @param posttx_exp scalar for average post-tx effects
#' @param assort_mix scalar (length 2 or length 9) for same/different class or all class relative mixing rates
#' @param long_traj logical for whether to include longitudinal risk trajectories in model
#' @param class_tx vector of class-specific relative treatment rates; if NULL, all classes have same treatment rates
#' @param beta_fit logical for whether model runs are being completed for beta fitting purposes; only age-specific prevalence
#' information for calculating likelihood will be output
#' @param summary_out logical for whether a smaller summary version of model output should be saved (TRUE) or if full model output
#' (all IDs/time steps) should be saved (FALSE); this will be overridden by beta_fit
#' 
init_model <- function(pop_dir = NULL,
                       total_pop = NULL,
                       out_id = NULL,
                       out_dir,
                       betas,
                       tx_rates,
                       mixing_df,
                       start_time,
                       end_time,
                       posttx_effects = TRUE,
                       posttx_exp = NULL,
                       assort_mix = NULL,
                       long_traj = TRUE,
                       class_tx = NULL,
                       beta_fit = FALSE,
                       summary_out = TRUE,
                       summary_month = FALSE,
                       age_brks = age_breaks
){
  
  if(is.null(pop_dir)){
    if(is.null(total_pop) | is.null(out_id)){
      print("must specify either total_pop or pop_dir to draw populations")
    }else{
      print("using total_pop to build populations")
      npop = length(out_id)
      pfile = FALSE
    }
  }else{
    print("using existing model output to run models")
    pfile=TRUE
    
    # check that files exist
    pop_files <- dir(pop_dir, full.names = TRUE)
    pop_num = as.numeric(gsub(".*_([0-9]+)*.parquet", "\\1", pop_files))
    npop = length(pop_num)
  }
  
  # # make sure `total_pop.parquet` and `total_pop_ids.Rdata` exist in pop_dir
  # if(!file.exists(paste0(pop_dir, "/total_pop_ids.Rdata"))){stop("total_pop_ids.Rdata does not exist")}
  # if(!file.exists(paste0(pop_dir, "/total_pop.parquet"))){stop("total_pop.parquet does not exist")}
  # total_pop <- read_parquet(paste0(pop_dir, "/total_pop.parquet"))
  # load(paste0(pop_dir, "/total_pop_ids.Rdata"))
  
  # make sure the output directory exists and is empty
  # if doesn't exist, create the directory
  if(dir.exists(out_dir)){
    if(length(dir(out_dir))>0) stop(glue("Output directory {out_dir} already containes files"))
  }else{
    dir.create(out_dir)
  }
  
  # if you're in the beta fitting module:
  if(beta_fit){
    # make sure if beta_fit=TRUE, summary_out=FALSE
    print("You specified beta fitting module. Summary output will not be saved")
    summary_out = FALSE
    
    # restrict out_id to first 25
    npop = 25
    
    # create a matrix to store the age-specific prevalences
    beta_out <- matrix(NA, nrow=10, ncol=npop)
  }
  
  # if there is class-assortative mixing, update the mixing matrix here
  if(!is.null(assort_mix)){
    
    # repeat the age categories for all class pairs
    tmp <- expand.grid(
      agec_cat = unique(mixing_df$agec_cat),
      agep_cat = unique(mixing_df$agep_cat),
      classp = 1:3,
      classc = 1:3
      )
    tmp <- tmp %>%
      # add back mij terms
      left_join(mixing_df, by=c("agec_cat", "agep_cat")) %>%
      # multiple by age_assort
      mutate(
        mij = ifelse(classp==classc, mij*assort_mix[1], mij*assort_mix[2]),
        classc = as.character(classc),
        classp = as.character(classp)
      ) %>%
      # create merging variables
      mutate(
        agec_cat_class = interaction(as.numeric(agec_cat), classc),
        agep_cat_class = interaction(as.numeric(agep_cat), classp)
      )
    mixing_df <- tmp
    rm(tmp)
  }
  
  
  # loop across each population (i.e., length of list of pop IDs)
  for(i in 1:npop){
    
    print(glue("iteration {i} start: {date()}"))
    
    # load the population
    if(pfile){
      pop <- read_parquet(pop_files[i])
      out_num = pop_num[i]
    }else{
      pop <- make_pop_from_ids(total_pop, out_id[[i]])
      out_num = i
    }

    # add betas and tx rates
    # these should be over-written no matter what
    pop <- pop[,which(!colnames(pop) %in% c("tx_rate", "beta"))]
    pop <- make_beta(pop, betas)
    pop <- make_tx_rates(pop, tx_rates)
    
    # if there is assortative mixing, add each individual's self-mixing
    # term based on age + class
    if(is.null(pop$indiv_mij)){
    if(!is.null(assort_mix)){
      pop <- pop %>%
        left_join(
          mixing_df %>%
          filter(agep_cat==agec_cat & classc==classp) %>%
          mutate(agec_cat = as.numeric(agec_cat), classc=as.numeric(classc)) %>%
          select(agec_cat, classc, mij) %>%
          rename(indiv_mij = mij),
        by=c("age_cat"="agec_cat", "class"="classc")
      )
    }else{
      # if there is no assortative mixing, just make self-mixing term
      # based on age categories
        pop <- pop %>%
          left_join(mixing_df %>% 
                      filter(agep_cat==agec_cat) %>% 
                      select(agec_cat, mij) %>%
                      rename(indiv_mij = mij), 
                    by=c("age_cat"="agec_cat"))
      }
    }
    
    # if long_traj=FALSE, remove effects of longitudinal trajectories
    if(!long_traj){
      pop$gamma = mean(pop$gamma, na.rm=T)
    }
    
    # if posttx_effects=FALSE, remove effects of post-tx reductions in risk
    if(!posttx_effects){
      pop <- pop %>%
        mutate(posttx_red = 0)
    }
    
    # if posttx_exp is non-null, adjust the post-tx reductions in risk
    if(!is.null(posttx_exp)){
      pop$posttx_red = pop$posttx_red^posttx_exp
    }
    
    # check whether class-specific treatment rates exist and, if so,
    # are they the right structure
    if(!is.null(class_tx) & length(class_tx)!=length(unique(pop$class))){
      stop("class_prob must be the same length as number of classes")
    }

    # run model
    out <- run_model(pop = pop, 
                     mixing_df = mixing_df, 
                     class_tx_prob = class_tx, 
                     assort_mix = assort_mix,
                     start_time = start_time, 
                     end_time = end_time)

    # create and store summary outputs from each simulation
    # by adding simulation number + creating indexed file name

    if(beta_fit){
      # make continued data-frame of age-specific prev
      beta_out[,i] <- make_age_prev(out, age_brks)
      if(i<4){
        make_out_summary(out, filename=glue("{out_dir}/summary_out_{out_num}.parquet"), sim_num=out_num)
      }
    } else if(summary_out){
      make_out_summary(out, filename=glue("{out_dir}/summary_out_{out_num}.parquet"), sim_num=out_num, monthly=summary_month)
      make_out_final_vis(out, filename=glue("{out_dir}/final_vis_out_{out_num}.parquet"), sim_num=out_num)
    }else{
      out$sim_num = out_num
      write_parquet(out, glue("{out_dir}/out_{out_num}.parquet"))
    }
    
    # clean-up before next population
    rm(out)
    rm(pop)
    
  }
  
  # if beta fitting, print out the data frame of all age-specific prevs
  if(beta_fit){
    colnames(beta_out) = glue("sim_num{1:npop}")
    write_parquet(as_tibble(beta_out), glue("{out_dir}/out_beta.parquet"))
  }
  
  # end
}


#' function to run transmission model for a given population set,
#' called from within init_model
#' @param pop df, population read in and modified in init_model
#' @param mixing_df df of age class i, j, and mixing intensity between ages, to be passed to make_infection function
#' @param start_time time point to start running the model; if >0, assumed that pop_dir points to prior model output
#' @param end_time time point to stop the model run, useful for calibrating betas (stop at t=9135)
#' @param class_tx_prob vector of class-specific relative treatment rates; if NULL, all classes have same treatment rates
run_model <- function(pop, 
                      mixing_df, 
                      class_tx_prob=NULL, 
                      assort_mix=NULL,
                      start_time, 
                      end_time){
  
  # outcome and state variables that will be used throughout model
  var_to_fill = c("oc_sero_pos", "oc_n_inf", "oc_n_new_inf", "oc_n_superinf", 
                  "oc_n_posttx_inf", "oc_n_posttx_new_inf", "oc_n_tx", "oc_n_tx_success", "oc_n_clearance", 
                  "oc_epsilon", "st_char", "t_since_infection")
  
  # if start_time > 0, the population object should have all outcome variables 
  # already created
  # stop model run if not
  if(start_time>0){
    stopifnot(
      var_to_fill %in% names(pop)
    )
  }
  
  # if start_time=0, add the following output variables:
  if(start_time==0){
    # sero_pos (ever infected)
    # n_inf (total number of infections during model period)
    # n_new_inf (number of infections moving from negative to acute infection stage)
    # n_superinf (number of infections moving from acute to acute infection stage or chronic to chronic infection stage)
    # n_posttx_inf (number of infections after treated at least once i.e. epsilon=1)
    # n_tx (number of treatment instances)
    # n_tx_success (number of successful treatment instances)
    # epsilon (indicator for whether ever treated)
    # n_clearance (number of times moving from acute to negative)
    #
    # and make values NA at all times after min(time) (entry time) for
    # each person, so previous time-point steps can be filled for
    # later time points
    entry_ind = which(pop$tau==0)
    
    pop$oc_sero_pos = NA_real_
    pop$oc_sero_pos[entry_ind] = pop$hcv_ab[entry_ind]

    pop$oc_n_inf = NA_real_
    pop$oc_n_new_inf = NA_real_
    pop$oc_n_superinf = NA_real_
    pop$oc_n_posttx_inf = NA_real_
    pop$oc_n_posttx_new_inf = NA_real_
    pop$oc_n_tx = NA_real_
    pop$oc_n_tx_success = NA_real_
    pop$oc_n_clearance = NA_real_
    pop$oc_epsilon = NA_real_
    
    pop$oc_n_inf[entry_ind] = 0
    pop$oc_n_new_inf[entry_ind] = 0
    pop$oc_n_superinf[entry_ind] = 0
    pop$oc_n_posttx_inf[entry_ind] = 0
    pop$oc_n_posttx_new_inf[entry_ind] = 0
    pop$oc_n_tx[entry_ind] = 0
    pop$oc_n_tx_success[entry_ind] = 0
    pop$oc_n_clearance[entry_ind] = 0
    pop$oc_epsilon[entry_ind] = 0

    # create and define infection states at entry into population
    # individuals can be in 
    # N - HCV negative, eligible for infection to move to st_acute
    # A - infected, eligible for clearance, superinfection, or moving to chronic infection
    # C - infected, eligible for treatment (successful or failure) or superinfection
    # at entry time, can either be N (if hcv_vir==0) or C (hcv_vir==1)
    # t_since_infection should be 0 for all, since no one is in acute state, NA at non-entry time points
    # again, states should be NA at all non-entry time points
    pop$st_char = NA_real_
    pop$st_char[pop$tau==0 & pop$hcv_vir==0] = 1
    pop$st_char[pop$tau==0 & pop$hcv_vir==1] = 3
    pop$t_since_infection = NA_real_
    pop$t_since_infection[entry_ind] = 0
  }
  
  # get model times to iterate across
  # remove those that are > end_time
  # and less than start_time
  times <- sort(unique(pop$time))
  times <- times[times<=end_time & times>=start_time]
  
  # make population into a matrix for model iterations/speed
  pop <- data.matrix(pop)
  
  # start at the second time point (i.e., states at first time point must already exist,
  # whether start_time=0 or some later time point)
  for(ii in 2:length(times)){
    
    # take the current time of interest + previous time point
    t_ind = times[ii]
    t_ind_prev = times[ii-1]
    
    # filter data to just time point of interest + previous time point
    # carry over the previous time points state + outcome variables
    st_prev = pop[pop[,"time"]==t_ind_prev,]
    
    # filter to just time point of interest
    # and those missing st_char (i.e., not entry at t_ind)
    tmp_pop_ind = which(pop[,"time"]==t_ind)
    tmp_pop_ind = tmp_pop_ind[which(is.na(pop[tmp_pop_ind,"st_char"]))]
    
    tmp <- pop[tmp_pop_ind,]
    
    # fill in state/outcome variables
    tmp[, var_to_fill] = st_prev[match(tmp[,"id"], st_prev[,"id"]), var_to_fill]
    tmp = cbind(tmp, st_char_prev=tmp[,"st_char"])
    tmp[,"st_char"] = NA_real_
    rm(st_prev)
    
    # ............................................
    # TREATMENT ..................................
    # ............................................
    
    tmp <- make_treatment(tmp, class_tx_prob)
    
    # ............................................
    # BECOME CHRONIC .............................
    # ............................................
    # only for those in acute class,
    # if t_since_infection > t_become_chronic, mark 
    # as moving into chronic stage
    tmp = cbind(tmp, tmp_chronic=0)
    tmp[tmp[,"t_since_infection"]>tmp[,"t_become_chronic"], "tmp_chronic"] = 1
    
    # ............................................
    # CLEARANCE ..................................
    # ............................................
    # only for those in the acute class, draw clearance event
    tmp = cbind(tmp, tmp_clearance = 0)
    ind = which(tmp[,"st_char_prev"]==2)
    tmp[ind,"tmp_clearance"] = rbinom(length(ind), 1, tmp[ind,"pr_clearance"])
    
    # ............................................
    # INFECTION ..................................
    # ............................................
    
    tmp <- make_infection(tmp, mixing_df, assort_mix)
    
    # ............................................
    # BOOK-KEEPING ...............................
    # ............................................
    
    # update individual states
    tmp <- update_state(tmp)
    
    # add tmp back into pop, replacing the times just run
    # and only those that don't already have st_char (i.e. entries into pop)
    pop[tmp_pop_ind, var_to_fill] = tmp[,var_to_fill]
   
    rm(tmp)
  }
  
  pop <- as.data.frame(pop)
  
  return(pop)
  
}



## ................................................
## OUTPUT FUNCTIONS............................####
## ................................................


#' function to create a summary data frame from model output,
#' including annual incidence, treatment rates, prevalence, etc.
#' This function is called if `summary_out=TRUE` in init_model
#' 
#' @param out df, output of model simulation
#' @param delta_t numeric, length of time step
#'

make_out_summary_int <- function(out, delta_t){
  
  # start with adding years
  out$year = floor(year0 + out$time/365.25)

  # prevalence in 2015, 2030 
  # seroprevalence + actively infected prevalence
  out_prev <- out %>%
    # take last visit for each year
    group_by(year, id) %>%
    slice(n()) %>%
    group_by(year) %>%
    # summarize sero prevalence + prevalence of active infection
    summarize(
      prev_sero = mean(oc_sero_pos),
      prev_inf = mean(st_char!=1),
      ever_tx = mean(oc_epsilon))
  
  # annual outcomes: 
  # incidence (all infections, new infections, person time = all; post-tx infections, person time = if oc_epsilon==1)
  # treatment (# doses; rate per PY, time in C; proportion treated at least once, denom=in C at least once)
  out_annual <- out %>% 
    arrange(time) %>% 
    group_by(id) %>%
    # take difference of oc_n variables to get incident values
    mutate(
      across(starts_with("oc_n"), function(x) c(0, diff(x))),
    ) %>% 
    group_by(year) %>%
    summarise(
      # get the sum of oc_n events
      across(starts_with("oc_n"), sum, na.rm=T),
      # get the follow-up time variables
      time_fu = delta_t*n(),
      time_in_C = delta_t*sum(st_char==3),
      time_in_A = delta_t*sum(st_char==2),
      time_in_N = delta_t*sum(st_char==1),
      time_posttx = delta_t*sum(oc_epsilon==1),
      time_posttx_new = delta_t*sum(oc_epsilon==1 & st_char==1),
      time_posttx_inf = delta_t*sum(oc_epsilon==1 & (st_char>1))
    ) %>%
    # calculate rates
    mutate(
      inc_oc_n_inf = oc_n_inf / time_fu * 365.25,
      inc_oc_n_new_inf = oc_n_new_inf / time_in_N * 365.25,
      inc_oc_n_superinf = oc_n_superinf / (time_in_C + time_in_A) * 365.25,
      inc_oc_n_posttx_inf = oc_n_posttx_inf / time_posttx * 365.25,
      inc_oc_n_posttx_new_inf = oc_n_posttx_new_inf / time_posttx_new * 365.25,
      inc_oc_n_tx = oc_n_tx / time_in_C * 365.25,
      inc_oc_n_tx_success = oc_n_tx_success / time_in_C * 365.25,
      inc_oc_n_clearance = oc_n_clearance / time_in_A * 365.25
    )
  
  
  # overall outcomes (all time):
  # incidence
  # treatment
  # average number of times treated per-person
  # average number of times infected per-person (overall, superinfection, post-tx)
  out_overall <- out_annual %>%
    ungroup %>%
    # sum the time/count outcomes across years
    summarize(across(starts_with("time") | starts_with("oc_n"), sum, na.rm=T)) %>%
    # calculate rates
    mutate(
      inc_oc_n_inf = oc_n_inf / time_fu * 365.25,
      inc_oc_n_new_inf = oc_n_new_inf / time_in_N * 365.25,
      inc_oc_n_superinf = oc_n_superinf / (time_in_C + time_in_A) * 365.25,
      inc_oc_n_posttx_inf = oc_n_posttx_inf / time_posttx * 365.25,
      inc_oc_n_posttx_new_inf = oc_n_posttx_new_inf / time_posttx_new * 365.25,
      inc_oc_n_tx = oc_n_tx / time_in_C * 365.25,
      inc_oc_n_tx_success = oc_n_tx_success / time_in_C * 365.25,
      inc_oc_n_clearance = oc_n_clearance / time_in_A * 365.25
    )
  
  # combine everything into one df for output
  out_overall$year = 9999
  out_summary <- left_join(out_prev, out_annual, by="year") 
  out_summary <- bind_rows(out_summary, out_overall)
  
  return(out_summary)
  
}


#' function to create a summary data frame from model output,
#' including MONTHLY incidence, treatment rates, prevalence, etc.
#' This function is called if `summary_out=TRUE` in init_model
#' 
#' @param out df, output of model simulation
#' @param delta_t numeric, length of time step
#'

make_out_summary_monthly_int <- function(out, delta_t){
  
  # start with adding years + months
  out$year = floor(year0 + out$time/365.25)
  out$month = ceiling((out$time - (out$year-year0)*365.25)/30.41667)
  out$month[out$month==0] = 1
  out$month[out$month==13] = 12
  
  # prevalence in 2015, 2030 
  # seroprevalence + actively infected prevalence
  out_prev <- out %>%
    # take last visit for each year
    group_by(month, year, id) %>%
    slice(n()) %>%
    group_by(month, year) %>%
    # summarize sero prevalence + prevalence of active infection
    summarize(
      prev_sero = mean(oc_sero_pos),
      prev_inf = mean(st_char!=1),
      ever_tx = mean(oc_epsilon))
  
  # annual outcomes: 
  # incidence (all infections, new infections, person time = all; post-tx infections, person time = if oc_epsilon==1)
  # treatment (# doses; rate per PY, time in C; proportion treated at least once, denom=in C at least once)
  out_annual <- out %>% 
    arrange(time) %>% 
    group_by(id) %>%
    # take difference of oc_n variables to get incident values
    mutate(
      across(starts_with("oc_n"), function(x) c(0, diff(x))),
    ) %>% 
    group_by(month, year) %>%
    summarise(
      # get the sum of oc_n events
      across(starts_with("oc_n"), sum, na.rm=T),
      # get the follow-up time variables
      time_fu = delta_t*n(),
      time_in_C = delta_t*sum(st_char==3),
      time_in_A = delta_t*sum(st_char==2),
      time_in_N = delta_t*sum(st_char==1),
      time_posttx = delta_t*sum(oc_epsilon==1),
      time_posttx_new = delta_t*sum(oc_epsilon==1 & st_char==1),
      time_posttx_inf = delta_t*sum(oc_epsilon==1 & (st_char>1))
    ) %>%
    # calculate rates
    mutate(
      inc_oc_n_inf = oc_n_inf / time_fu * 365.25,
      inc_oc_n_new_inf = oc_n_new_inf / time_in_N * 365.25,
      inc_oc_n_superinf = oc_n_superinf / (time_in_C + time_in_A) * 365.25,
      inc_oc_n_posttx_inf = oc_n_posttx_inf / time_posttx * 365.25,
      inc_oc_n_posttx_new_inf = oc_n_posttx_new_inf / time_posttx_new * 365.25,
      inc_oc_n_tx = oc_n_tx / time_in_C * 365.25,
      inc_oc_n_tx_success = oc_n_tx_success / time_in_C * 365.25,
      inc_oc_n_clearance = oc_n_clearance / time_in_A * 365.25
    )
  
  
  # overall outcomes (all time):
  # incidence
  # treatment
  # average number of times treated per-person
  # average number of times infected per-person (overall, superinfection, post-tx)
  out_overall <- out_annual %>%
    ungroup %>%
    # sum the time/count outcomes across years
    summarize(across(starts_with("time") | starts_with("oc_n"), sum, na.rm=T)) %>%
    # calculate rates
    mutate(
      inc_oc_n_inf = oc_n_inf / time_fu * 365.25,
      inc_oc_n_new_inf = oc_n_new_inf / time_in_N * 365.25,
      inc_oc_n_superinf = oc_n_superinf / (time_in_C + time_in_A) * 365.25,
      inc_oc_n_posttx_inf = oc_n_posttx_inf / time_posttx * 365.25,
      inc_oc_n_posttx_new_inf = oc_n_posttx_new_inf / time_posttx_new * 365.25,
      inc_oc_n_tx = oc_n_tx / time_in_C * 365.25,
      inc_oc_n_tx_success = oc_n_tx_success / time_in_C * 365.25,
      inc_oc_n_clearance = oc_n_clearance / time_in_A * 365.25
    )
  
  # combine everything into one df for output
  out_overall$year = 9999
  out_overall$month = 9999
  out_summary <- left_join(out_prev, out_annual, by=c("month", "year"))
  out_summary <- bind_rows(out_summary, out_overall)
  
  return(out_summary)
  
}

#' function to make output summary for each class
#' and overall, using make_out_summary
#' Directly write file
#' 
#' @param out df, output of model simulation
#' @param filename, path to write output
#' @param sim_num, numeric, which population/simulation is being output
#' @param delta_t numeric, length of time step
#'

make_out_summary <- function(out, filename, sim_num, delta_t=dt, monthly=FALSE){
  
  if(monthly){
    tmp_overall = make_out_summary_monthly_int(out, delta_t)
    tmp_overall$class = "All"
    
    # for each class
    tmp_class1 = make_out_summary_monthly_int(out[out$class==1,], delta_t)
    tmp_class2 = make_out_summary_monthly_int(out[out$class==2,], delta_t)
    tmp_class3 = make_out_summary_monthly_int(out[out$class==3,], delta_t)    
  }else{
    # make overall
    tmp_overall = make_out_summary_int(out, delta_t)
    tmp_overall$class = "All"
    
    # for each class
    tmp_class1 = make_out_summary_int(out[out$class==1,], delta_t)
    tmp_class2 = make_out_summary_int(out[out$class==2,], delta_t)
    tmp_class3 = make_out_summary_int(out[out$class==3,], delta_t)
  }
  
  tmp_class1$class = "1"
  tmp_class2$class = "2"
  tmp_class3$class = "3"
  
  out <- bind_rows(tmp_overall, tmp_class1, tmp_class2, tmp_class3)
  
  out$sim_num = sim_num
  
  write_parquet(out, filename)
  
}

#' function to store the last visit for each individual
#' 
#' @param out df, output of model simulation
#' @param filename, path to write output
#' @param sim_num, numeric, which population/simulation is being output
#'

make_out_final_vis <- function(out, filename, sim_num){
  
  # make overall  
  tmp <- out %>%
    group_by(id) %>%
    slice(n())
  tmp$sim_num = sim_num
  
  write_parquet(tmp, filename)
  
}

#' function to create vector of age-class specific
#' prevalence for use in beta fitting module
#' 
#' @param out df, output of model simulation
#' 

make_age_prev <- function(out, age_breaks){
  
  # start with adding years
  # out$year_num = year0 + out$time/365.25
  # out$year = floor(out$year_num)
  
  # filter out the last visit for each ID in 2014
  # and calculate age specific seroprevalence
  # tmp <- out[out$year==2015,]
  # tmp <- tmp %>% group_by(id) %>% slice(n())
  
  # filter out the end of 2014 visit (fixed time point for beta estimation)
  # t=9135
  tmp <- out[out$time==9135,]
  tmp$age_cat <- cut(tmp$age, age_breaks, right=FALSE)
  out <- tapply(tmp$oc_sero_pos, tmp$age_cat, mean)
  
  return(out)
}



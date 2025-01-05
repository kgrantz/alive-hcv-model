# Model of HCV transmission among PWID

Individual-based model of HCV transmission among people who inject drugs, as presented in Grantz et al., "Impact of risk heterogeneity on the feasibility of Hepatitis C elimination among people who inject drugs: a modeling study".

There are two general steps to running this model:

1. Creating the population of people who inject drugs (`create-population.R`)
2. Running HCV transmission simulations (`run-model-example`)

Output from simulations can be used for parameter calibration using grid searches and/or forward simulations of specific transmission and treatment scenarios.

The provided parameter file includes examples of the necessary parameters, but values have been changed to preserve confidentiality of ALIVE participants. Original parameters can be provided upon request to repository owner.

Required parameters include:

- `age_enrol`: a named list of ages of enrollment for each risk class, in years
- `age_idu`: a named list of ages of initiation of injection drug use for each risk class, in years
- `bsero_byage`: df of age-specific seroprevalence
- `classes`: list with items `levels`, listing the labels used for each risk class, and `probs`, the probability of belonging to a certain risk class
- `mij`: long-form df of mixing matrix with relative rates of contact `mij` between contacts of age `agec_cat` and `agep_cat`
- `mspl_rf`: named list with data frame `times`, containing the days since enrollment at which the risk multiplier was evaluated, and matrix `pred`, containing the predicted median and 95% confidence intervals of the risk multiplier for each risk class.
- `mu`: df with all-cause mortality rates by `age_cat` and `year_cat`.

- `dt`, `nt`: time step and maximum time for running model, in days
- `hcv_age_param`: slope of the exponential curve of seropositivity by age
- `hcv_prev_scaled`: class-specific scalar on HCV prevalence, used in combination with `hcv_age_param`
- `hcv_vir_prob`: probability of being viremic if HCV seropositive at enrollment
- `n0`: starting population size; for testing, set to 100
- `nsim`: number of stochastic iterations to run of model; for testing, set to 5
- `sc_alpha`, `sc_beta`: parameters of beta distribution of probability of spontaneous clearance
- `sctime_min`, `sctime_max`: parameters of uniform distirbution of time from infection to spontaneous clearance
- `tx_eff`: treatment efficacy, or probability of clearance if treated
- `txred_alpha`, `txred_beta`: parameters of beta distribution of post-treatment reduction in infection risk
- various helper values and vectors for time-keeping (`age_breaks`, `age_cat`, `ages`, `year_breaks`, `year_labels`, `year0`)

#------------------------------------------------------------------------------#
# Replication of Brookhart M.A. et al. (2006)
# Variable Selection for Propensity Score Models.
# American Journal of Epidemiology, 163(12), 1149–1156.
#
# Replicator: K Luijken
# Co-pilot: B B L Penning de Vries
#
# Helper functions to run analyses
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Seeds ----
#------------------------------------------------------------------------------#
# Store seeds for each repetition for experiment 1, sensitivity analysis of
# experiment 1 and experiment 2

get_seeds <- function(rep){
  n_seed <- nrow(main_scenarios_1) * rep
  set.seed(08131)
  seeds_main_1 <- sample(1:1e8,
                         size = n_seed,
                         replace = FALSE)

  n_seed <- nrow(sensitivity_1) * rep
  set.seed(08132)
  seeds_sensitivity_1 <- sample(1:1e8,
                                size = n_seed,
                                replace = FALSE)

  n_seed <- nrow(main_scenarios_2) * rep
  set.seed(08133)
  seeds_main_2 <- sample(1:1e8,
                         size = n_seed,
                         replace = FALSE)

  return(list(seeds_main_1 = seeds_main_1, seeds_sensitivity_1 = seeds_sensitivity_1, seeds_main_2 = seeds_main_2))
}


#------------------------------------------------------------------------------#
# Single repetition helpers ----
#------------------------------------------------------------------------------#

# Experiment 1 ----
experiment1_single_rep <- function(datagen_scenario, ps_model_covs, seed){
  # Generate data
  data <- gen_data(seed = seed,
                   nobs = datagen_scenario[['nobs']],
                   sd_X1 = datagen_scenario[['sd_X1']],
                   sd_X2 = datagen_scenario[['sd_X2']],
                   sd_X3 = datagen_scenario[['sd_X3']],
                   beta0 = datagen_scenario[['beta0']],
                   beta1 = datagen_scenario[['beta1']],
                   beta2 = datagen_scenario[['beta2']],
                   beta3 = datagen_scenario[['beta3']],
                   alpha0 = datagen_scenario[['alpha0']],
                   alpha1 = datagen_scenario[['alpha1']],
                   alpha2 = datagen_scenario[['alpha2']],
                   alpha3 = datagen_scenario[['alpha3']],
                   alpha4 = datagen_scenario[['alpha4']],
                   trans = function(x) 1/(1+exp(-3 * x)) - 0.5)

  # Estimate crude effect
  effect_crude     <- glm(Y ~ A, data = data, family = poisson)$coefficients["A"]

  # Initialize results data.table
  results <- data.table::data.table(seed = numeric(),
                                    ps_model = character(),
                                    effect_crude = numeric(),
                                    effect_splines = numeric(),
                                    effect_quintiles = numeric(),
                                    c_stat = numeric())

  warnings <- data.table::data.table(ps_model = character(),
                         ps_warning = character(),
                         splines_warning = character(),
                         quintiles_1_warning = character(),
                         quintiles_2_warning = character(),
                         quintiles_3_warning = character(),
                         quintiles_4_warning = character(),
                         quintiles_5_warning = character())

  # Loop over different propensity score specifications
  for(i in 1:length(ps_model_covs)){
    # Estimate propensity score
    prop_score <- estimate_ps(ps_model_covs = ps_model_covs[[i]],
                      data = data)

    # Estimate outcome model
    effect_splines   <- estimate_effect_spline(PS = prop_score$PS,
                                               data = data)

    effect_quintiles <- estimate_effect_quintiles(PS = prop_score$PS,
                                                  data = data)

    # Results, saved by ps_model
    results <- rbind(results, data.table::data.table(seed = seed,
                                                     ps_model = paste(ps_model_covs[[i]], collapse = "+"),
                                                     effect_crude = effect_crude,
                                                     effect_splines = effect_splines$effect_A,
                                                     effect_quintiles = effect_quintiles$effect_A,
                                                     c_stat = prop_score$cstat_ps))
    warnings <- rbind(warnings, data.table::data.table(ps_model = paste(ps_model_covs[[i]], collapse = "+"),
                                                       ps_warning = as.character(prop_score$warning),
                                                       splines_warning = as.character(effect_splines$warning),
                                                       quintiles_1_warning = as.character(effect_quintiles$warning[1]),
                                                       quintiles_2_warning = as.character(effect_quintiles$warning[2]),
                                                       quintiles_3_warning = as.character(effect_quintiles$warning[3]),
                                                       quintiles_4_warning = as.character(effect_quintiles$warning[4]),
                                                       quintiles_5_warning = as.character(effect_quintiles$warning[5])))
  }

  return(list(results = results, warnings = warnings))
}

# Experiment 2 ----
experiment2_single_rep <- function(datagen_scenario,
                                   ps_model_covs,
                                   seed){
  # Generate data
  data <- gen_data(seed = seed,
                   nobs = datagen_scenario[['nobs']],
                   sd_X1 = datagen_scenario[['sd_X1']],
                   sd_X2 = datagen_scenario[['sd_X2']],
                   sd_X3 = datagen_scenario[['sd_X3']],
                   beta0 = datagen_scenario[['beta0']],
                   beta1 = datagen_scenario[['beta1']],
                   beta2 = datagen_scenario[['beta2']],
                   beta3 = datagen_scenario[['beta3']],
                   alpha0 = datagen_scenario[['alpha0']],
                   alpha1 = datagen_scenario[['alpha1']],
                   alpha2 = datagen_scenario[['alpha2']],
                   alpha3 = datagen_scenario[['alpha3']],
                   alpha4 = datagen_scenario[['alpha4']],
                   trans = identity)

  # Estimate propensity score
  prop_score <- estimate_ps(ps_model_covs = ps_model_covs,
                      data = data)

  # Estimate outcome model
  effect_crude     <- glm(Y ~ A, data = data, family = poisson)$coefficients["A"]
  effect_splines   <- estimate_effect_spline(PS = prop_score$PS,
                                             data = data)
  effect_quintiles <- estimate_effect_quintiles(PS = prop_score$PS,
                                                data = data)

  results <- data.table::data.table(effect_crude = effect_crude,
                                    effect_splines = effect_splines$effect_A,
                                    c_stat = prop_score$cstat_ps,
                                    effect_quintiles = effect_quintiles$effect_A)

  warnings <- data.table::data.table(ps_warning = as.character(prop_score$warning),
                                     splines_warning = as.character(effect_splines$warning),
                                     quintiles_1_warning = as.character(effect_quintiles$warning[1]),
                                     quintiles_2_warning = as.character(effect_quintiles$warning[2]),
                                     quintiles_3_warning = as.character(effect_quintiles$warning[3]),
                                     quintiles_4_warning = as.character(effect_quintiles$warning[4]),
                                     quintiles_5_warning = as.character(effect_quintiles$warning[5]))


  return(list(results = results, warnings = warnings))


}


#------------------------------------------------------------------------------#
# Workhorse ----
#------------------------------------------------------------------------------#
# Function to loop over repetitions within scenario

run_scenario <- function(use_datagen_scenario, ps_model_covs, seeds, rep, filepath,function_experiment){

  output <- NULL
  warnings <- NULL
  for(i in 1:rep){
    run <- function_experiment(datagen_scenario = use_datagen_scenario,
                                      ps_model_covs = ps_model_covs,
                                      seed = seeds[rep * (use_datagen_scenario[['scen_num']]-1) + i])

    output <- rbind(output, run$results)
    warnings <- rbind(warnings, run$warnings)
  }

  saveRDS(output, file = paste0(filepath,"S",use_datagen_scenario[['scen_num']],".rds"))
  saveRDS(warnings, file = paste0(filepath,"S",use_datagen_scenario[['scen_num']],"_warnings.rds"))
}

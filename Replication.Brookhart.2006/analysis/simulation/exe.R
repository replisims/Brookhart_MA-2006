#------------------------------------------------------------------------------#
# Replication of Brookhart M.A. et al. (2006)
# Variable Selection for Propensity Score Models.
# American Journal of Epidemiology, 163(12), 1149–1156.
#
# Replicator: K Luijken
# Co-pilot: B B L Penning de Vries
#
# Script to execute the simulation study
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Libraries and source code ----
#------------------------------------------------------------------------------#
# Load all functions to replicate the simulation (to be found in ./R folder)
devtools::load_all()

# Source helper scripts
source("./analysis/simulation/scenarios.R")
source("./analysis/simulation/run_sim.R")
source("./analysis/simulation/tables.R")
source("./analysis/simulation/graphs.R")

#------------------------------------------------------------------------------#
# Execute simstudy ----
#------------------------------------------------------------------------------#

# Main analysis experiment 1
rep <- 1000
all_seeds <- get_seeds(rep = rep)

apply(main_scenarios_1,
      MARGIN = 1,
      FUN = run_scenario,
      ps_model_covs = ps_covs_1,
      rep = rep,
      seeds = all_seeds$seeds_main_1,
      function_experiment = experiment1_single_rep,
      filepath = "./analysis/data/raw_data/experiment_1_main/")

apply(sensitivity_1,
      MARGIN = 1,
      FUN = run_scenario,
      ps_model_covs = ps_covs_1,
      rep = rep,
      seeds = all_seeds$seeds_sensitivity_1,
      function_experiment = experiment1_single_rep,
      filepath = "./analysis/data/raw_data/experiment_1_sensitivity/")

apply(main_scenarios_2,
      MARGIN = 1,
      FUN = run_scenario,
      ps_model_covs = "X1",
      rep = rep,
      seeds = all_seeds$seeds_main_2,
      function_experiment = experiment2_single_rep,
      filepath = "./analysis/data/raw_data/experiment_2_main/")


#------------------------------------------------------------------------------#
# Summarize results into tables ----
#------------------------------------------------------------------------------#
# Make sure to use correct values for important parameters
alpha4 <- 0.5

# Functions that produce tables as in original manuscript (latex format)

create_table1()
create_table2()
create_table3()

#------------------------------------------------------------------------------#
# Visualisations ----
#------------------------------------------------------------------------------#
# Make sure key parameters are set to correct values

rep <- 1000
alpha1 <- exp(seq(0, 0.2, by = 0.01)) # exponentiate values to obtain RRs
alpha4 <- 0.5
beta1  <- exp(seq(0, 1.25, by = 0.05))

# Functions that produce figures similar to original manuscript. Figures are
# saved as png files in ./analysis/figures/

create_figure2()
create_figure3()
create_figure4()

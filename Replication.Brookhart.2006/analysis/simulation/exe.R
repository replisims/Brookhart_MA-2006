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

# Set repetitions and seeds
rep <- 1000
all_seeds <- get_seeds(rep = rep)


# Main analysis experiment 1
apply(scenarios_1_main,
      MARGIN = 1,
      FUN = run_scenario,
      ps_model_covs = ps_covs_1,
      rep = rep,
      seeds = all_seeds$seeds_main_1,
      function_experiment = experiment1_single_rep,
      filepath = "./analysis/data/raw_data/experiment_1_main/")

# Sensitivity analysis experiment 1
apply(sensitivity_1,
      MARGIN = 1,
      FUN = run_scenario,
      ps_model_covs = ps_covs_1,
      rep = rep,
      seeds = all_seeds$seeds_sensitivity_1,
      function_experiment = experiment1_single_rep,
      filepath = "./analysis/data/raw_data/experiment_1_sensitivity/")

# Main analysis experiment 2
apply(scenarios_2_main,
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
# Functions that produce tables as in original manuscript (latex format)
alpha4 <- 0.5 # Make sure to use correct values for important parameters
create_table1()
create_table2()
create_table3()

#------------------------------------------------------------------------------#
# Visualisations ----
#------------------------------------------------------------------------------#
# Make sure key parameters are set to correct values
alpha4 <- 0.5

# Functions that produce figures similar to original manuscript. Figures are
# saved as png files in ./analysis/figures/

create_figure2(output = scenarios_2_main)
create_figure3(output = scenarios_2_main,
               alpha4 = alpha4)
create_figure4(output = scenarios_2_main,
               alpha4 = alpha4)

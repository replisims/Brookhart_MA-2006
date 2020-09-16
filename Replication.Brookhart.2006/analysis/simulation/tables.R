#------------------------------------------------------------------------------#
# Replication of Brookhart M.A. et al. (2006)
# Variable Selection for Propensity Score Models.
# American Journal of Epidemiology, 163(12), 1149–1156.
#
# Replicator: K Luijken
# Co-pilot: B B L Penning de Vries
#
# Script to generate output tables
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Helper functions ----
#------------------------------------------------------------------------------#
bias     <- function(x) mean(x - alpha4)
variance <- function(x) (1/(rep-1)) * sum((x - mean(x))^2)
mse      <- function(x) mean((x - alpha4)^2)


single_scenario <- function(filepath,
                            output_effect,
                            ps_covs_1 = ps_covs_1,
                            multiplication){
  results <- readRDS(paste0("./analysis/data/raw_data/", filepath))
  results <- data.frame(results)

  # Bias
  output_bias   <- rbind(by(results[,output_effect], results[,"ps_model"], function(x) bias(x) * multiplication))
  output_bias   <- cbind(output_bias, bias(results[,"effect_crude"]) * multiplication)

  # Variance
  output_var    <- rbind(by(results[,output_effect], results[,"ps_model"], function(x) variance(x) * multiplication))
  output_var    <- cbind(output_var, variance(results[,"effect_crude"]) * multiplication)

  # MSE
  output_mse    <- rbind(by(results[,output_effect], results[,"ps_model"], function(x) mse(x) * multiplication))
  output_mse    <- cbind(output_mse, mse(results[,"effect_crude"]) * multiplication)

  output        <- rbind(output_bias, output_var, output_mse)

  return(output)
}


add_c_stat <- function(filepath){
  results <- readRDS(paste0("./analysis/data/raw_data/", filepath))
  results <- data.frame(results)

  c_stat <- rbind(by(results[,"c_stat"], results[,"ps_model"], function(x) mean(x)))
  return(c_stat)
}

#------------------------------------------------------------------------------#
# Table 1: Experiment 1, spline
#------------------------------------------------------------------------------#

create_table1 <- function(){
  table1 <- xtable::xtable(rbind(rep(NA, times = 8),
                        single_scenario(filepath = "experiment_1_main/S1.rds",
                                        output_effect = "effect_splines",
                                        ps_covs_1 = ps_covs_1,
                                        multiplication = 10),
                        c(add_c_stat(filepath = "experiment_1_main/S1.rds"),NA),
                        rep(NA, times = 8),
                        single_scenario(filepath = "experiment_1_main/S2.rds",
                                        output_effect = "effect_splines",
                                        ps_covs_1 = ps_covs_1,
                                        multiplication = 100),
                        c(add_c_stat(filepath = "experiment_1_main/S2.rds"),NA)),digits = 2)
  rownames(table1) <- c("n = 500","Bias x 10","Variance x 10","MSE x 10","Average c-statistic",
                        "n = 2,500","Bias x 100","Variance x 100","MSE x 100","Average cstatistic")

  addtorow <- list()
  addtorow$pos <- list(0,0)
  addtorow$command <- c("& \\multicolumn{8}{c}{Variable(s) in propensity score model} \\\\\n",
                        "& X1 & X1+X2 & X1+X2+X3 & X1+X3 & X2 & X2+X3 & X3 & None \\\\\n")

  print(table1, add.to.row = addtorow, include.colnames = FALSE)
}
  #------------------------------------------------------------------------------#
# Table 2: Experiment 1, subclassification
#------------------------------------------------------------------------------#

create_table2 <- function(){
  table2 <- xtable::xtable(rbind(rep(NA, times = 8),
                         single_scenario(filepath = "experiment_1_main/S1.rds",
                                         output_effect = "effect_quintiles",
                                         ps_covs_1 = ps_covs_1,
                                         multiplication = 10),
                         rep(NA, times = 8),
                         single_scenario(filepath = "experiment_1_main/S2.rds",
                                         output_effect = "effect_quintiles",
                                         ps_covs_1 = ps_covs_1,
                                         multiplication = 100)),digits = 2)
  rownames(table2) <- c("n = 500","Bias x 10","Variance x 10","MSE x 10",
                        "n = 2,500","Bias x 100","Variance x 100","MSE x 100")

  addtorow <- list()
  addtorow$pos <- list(0,0)
  addtorow$command <- c("& \\multicolumn{8}{c}{Variable(s) in propensity score model} \\\\\n",
                        "& X1 & X1+X2 & X1+X2+X3 & X1+X3 & X2 & X2+X3 & X3 & None \\\\\n")

  print(table2, add.to.row = addtorow, include.colnames = FALSE)
}

#------------------------------------------------------------------------------#
# Table 3: Experiment 1, sensitivity
#------------------------------------------------------------------------------#

create_table3 <- function(){
  table3 <- rbind(rep(NA, times = 8),
                  single_scenario(filepath = "experiment_1_main/S1.rds",
                                  output_effect = "effect_splines",
                                  ps_covs_1 = ps_covs_1,
                                  multiplication = 10))

  for(i in 1:nrow(sensitivity_1)){
    table3 <- rbind(table3,
                    rep(NA, times = 8),
                    single_scenario(filepath = paste0("experiment_1_sensitivity/S",i,".rds"),
                                    output_effect = "effect_splines",
                                    ps_covs_1 = ps_covs_1,
                                    multiplication = 10))
  }
  table3 <- round(table3, digits = 2)

  table3 <- cbind(c("Original","Bias x 10","Variance x 10","MSE x 10",
                    "Decrease in the variance of X1","Bias x 10","Variance x 10","MSE x 10",
                    "Increase in the variance of X1","Bias x 10","Variance x 10","MSE x 10",
                    "Decrease in the variance of X2","Bias x 10","Variance x 10","MSE x 10",
                    "Increase in the variance of X2","Bias x 10","Variance x 10","MSE x 10",
                    "Decrease in the variance of X3","Bias x 10","Variance x 10","MSE x 10",
                    "Increase in the variance of X3","Bias x 10","Variance x 10","MSE x 10",
                    "Decrease in alpha4","Bias x 10","Variance x 10","MSE x 10",
                    "Increase in alpha4","Bias x 10","Variance x 10","MSE x 10",
                    "Decrease in beta0","Bias x 10","Variance x 10","MSE x 10"), table3)

  table3 <- xtable::xtable(table3,digits = 2)

  addtorow <- list()
  addtorow$pos <- list(0,0)
  addtorow$command <- c("& \\multicolumn{8}{c}{Variable(s) in propensity score model} \\\\\n",
                        "& X1 & X1+X2 & X1+X2+X3 & X1+X3 & X2 & X2+X3 & X3 & None \\\\\n")

  print(table3, add.to.row = addtorow, include.colnames = FALSE)
}


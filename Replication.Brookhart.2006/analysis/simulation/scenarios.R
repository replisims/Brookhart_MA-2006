#------------------------------------------------------------------------------#
# Replication of Brookhart M.A. et al. (2006)
# Variable Selection for Propensity Score Models.
# American Journal of Epidemiology, 163(12), 1149–1156.
#
# Replicator: K Luijken
# Co-pilot: B B L Penning de Vries
#
# Script to generate scenarios
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Experiment 1 ----
#------------------------------------------------------------------------------#

# Simulation parameter values
nobs   <- c(500, 2500)
sd_X1  <- 1
sd_X2  <- 1
sd_X3  <- 1
alpha0 <- 0.5
alpha1 <- 4
alpha2 <- 1
alpha3 <- 0
alpha4 <- 0.5
beta0  <- 0
beta1  <- 0.5
beta2  <- 0
beta3  <- 0.75

# Scenario grid main analysis
main_scenarios_1 <- expand.grid(nobs = nobs,
                                sd_X1 = sd_X1,
                                sd_X2 = sd_X2,
                                sd_X3 = sd_X3,
                                alpha0 = alpha0,
                                alpha1 = alpha1,
                                alpha2 = alpha2,
                                alpha3 = alpha3,
                                alpha4 = alpha4,
                                beta0 = beta0,
                                beta1 = beta1,
                                beta2 = beta2,
                                beta3 = beta3)

main_scenarios_1$scen_num <- seq(1, nrow(main_scenarios_1), by = 1)

# Scenario grid sensitivity analyses
default <- data.frame(nobs = 500, sd_X1 = sd_X1, sd_X2 = sd_X2, sd_X3 = sd_X3,
                      alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, alpha4 = alpha4,
                      beta0 = beta0, beta1 = beta1, beta2 = beta2, beta3 = beta3)
## Set defaul values
sensitivity_1 <- default[rep(seq_len(nrow(default)), each=9),]

## Change one parameter at a time
sensitivity_1$sd_X1[1] <- 0.5
sensitivity_1$sd_X1[2] <- 1.5
sensitivity_1$sd_X2[3] <- 0.5
sensitivity_1$sd_X2[4] <- 1.5
sensitivity_1$sd_X3[5] <- 0.5
sensitivity_1$sd_X3[6] <- 1.5
sensitivity_1$alpha4[7] <- 0.25
sensitivity_1$alpha4[8] <- 1
sensitivity_1$beta0[9] <- -1
sensitivity_1$scen_num <- seq(1, nrow(sensitivity_1), by = 1)

# Propensity models main analysis
ps_covs_1 <- list(model1 = "X1",
                    model2 = "X2",
                    model3 = "X3",
                    model4 = c("X1", "X2"),
                    model5 = c("X1", "X3"),
                    model6 = c("X2", "X3"),
                    model7 = c("X1", "X2", "X3"))

#------------------------------------------------------------------------------#
# Experiment 2 ----
#------------------------------------------------------------------------------#

# Simulation parameter values
nobs   <- c(500, 2500)
sd_X1  <- 1
alpha0 <- 0.5 # Not reported in manuscript
alpha1 <- seq(0, 0.2, by = 0.01)
alpha4 <- 0.5
beta0  <- 0   # Not reported in manuscript
beta1  <- seq(0, 1.25, by = 0.05)

# Scenario grid main analysis
main_scenarios_2 <- expand.grid(nobs = nobs,
                                sd_X1 = sd_X1,
                                sd_X2 = 0,
                                sd_X3 = 0,
                                alpha0 = alpha0,
                                alpha1 = alpha1,
                                alpha2 = 0,
                                alpha3 = 0,
                                alpha4 = alpha4,
                                beta0 = beta0,
                                beta1 = beta1,
                                beta2 = 0,
                                beta3 = 0)
main_scenarios_2$scen_num <- seq(1, nrow(main_scenarios_2), by = 1)

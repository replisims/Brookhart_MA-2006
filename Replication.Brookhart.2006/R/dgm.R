#------------------------------------------------------------------------------#
# Replication of Brookhart M.A. et al. (2006)
# Variable Selection for Propensity Score Models.
# American Journal of Epidemiology, 163(12), 1149–1156.
#
# Replicator: K Luijken
# Co-pilot: B B L Penning de Vries
#
# Script to generate data
#------------------------------------------------------------------------------#


#' Generate simulation data for replication Brookhart et al. (2006)
#'
#' @param seed The seed used to generate a dataset.
#' @param nobs The number of observations of the generated dataset.
#' @param sd_X1 The standard deviation of covariate X1.
#' @param sd_X2 The standard deviation of covariate X2.
#' @param sd_X3 The standard deviation of covariate X3.
#' @param beta0 The intercept of the model for the conditional mean of exposure A on a log-scale.
#' @param beta1 The conditional association between covariate X1 and exposure A on a log-scale.
#' @param beta2 The conditional association between covariate X2 and exposure A on a log-scale.
#' @param beta3 The conditional association between covariate X3 and exposure A on a log-scale.
#' @param alpha0 The intercept of the model for the conditional mean of outcome Y on a log-scale.
#' @param alpha1 The conditional association between covariate X1 and outcome Y on a log-scale.
#' @param alpha2 The conditional association between covariate X2 and outcome Y on a log-scale.
#' @param alpha3 The conditional association between covariate X3 and outcome Y on a log-scale.
#' @param alpha4 The conditional association between exposure A and outcome Y on a log-scale. Set to 0.5 for all scenarios in this replication.
#' @param trans Transformations of covariate X1 in generating the outcome, Y. For experiment 1, use trans = identity; for experiment 2, use trans = function(x) 1/(1+exp(-3 * x)) * -.5
#' on a log-scale.
#'
#' @return The function returns a simulated dataset containing a dichotomous exposure, an outcome with a Poisson distribution, and continuous independent confounders. Both simulation experiment 1 and experiment 2 of the manuscript employ the same data-generating process.
gen_data <- function(seed,
                     nobs,
                     sd_X1,
                     sd_X2,
                     sd_X3,
                     beta0,
                     beta1,
                     beta2,
                     beta3,
                     alpha0,
                     alpha1,
                     alpha2,
                     alpha3,
                     alpha4 = 0.5,
                     trans = identity){
  # Generate data
  set.seed(seed)

  # Generate covariates
  X1 <- rnorm(nobs, mean = 0, sd = sd_X1)
  X2 <- rnorm(nobs, mean = 0, sd = sd_X2)
  X3 <- rnorm(nobs, mean = 0, sd = sd_X3)

  # Genarate exposure A
  cond_mean_A <- pnorm(beta0 + X1 * beta1 + X2 * beta2 + X3 * beta3)
  A  <- rbinom(nobs, 1, cond_mean_A)

  # Generate outcome Y
  cond_mean_Y <- exp(alpha0 + trans(X1) * alpha1 + X2 * alpha2 +
                       X3 * alpha3 + A * alpha4)
  Y  <- rpois(nobs, cond_mean_Y)

  # Store as dataframe
  out_df <- data.frame(Y = Y, A = A, X1 = X1, X2 = X2, X3 = X3)

  return(out_df)
}


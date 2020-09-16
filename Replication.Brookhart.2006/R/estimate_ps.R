#------------------------------------------------------------------------------#
# Replication of Brookhart M.A. et al. (2006)
# Variable Selection for Propensity Score Models.
# American Journal of Epidemiology, 163(12), 1149–1156.
#
# Replicator: K Luijken
# Co-pilot: B B L Penning de Vries
#
# Helper function to estimate propensity score
#------------------------------------------------------------------------------#


#' Estimating the propensity score model for replication Brookhart et al. (2006)
#'
#' @param ps_model_covs Character vector of covariates that should be included in the propensity score model.
#' @param data The dataset containing the exposure variable A, outcome variable Y and covariates.
#'
#' @return Provides a propensity score and c-statistic of the propensity score model. Stores warnings that occurred during propensity model fitting.
#'
#' @import Hmisc
#' @import simsalapar
estimate_ps <- function(ps_model_covs, data){
  PS_formula <- as.formula(paste0("A ~ ", paste(ps_model_covs, collapse = "+")))

  # Estimate propensity score model and store possible warnings using simsalapar::tryCatch.W.E
  PS_mod     <- simsalapar::tryCatch.W.E(
                    glm(PS_formula, data = data, family = binomial(link = "probit")))

  # Estimate propensity scores
  PS         <- predict(object = PS_mod$value,
                        newdata = data[,c('A',ps_model_covs)],
                        type = "response")

  # Store c-statistic (estimated using Hmisc package)
  cstat_ps   <- Hmisc::somers2(PS, data[,"A"])["C"]

  # When no warning occured, save NA in warning output matrix
  if(is.null(PS_mod$warning)){PS_mod$warning <- NA}

  return(list(PS = PS, cstat_ps = cstat_ps, warning = PS_mod$warning))
}

#------------------------------------------------------------------------------#
# Replication of Brookhart M.A. et al. (2006)
# Variable Selection for Propensity Score Models.
# American Journal of Epidemiology, 163(12), 1149–1156.
#
# Replicator: K Luijken
# Co-pilot: B B L Penning de Vries
#
# Helper function to estimate outcome model including the ps as a spline function
#------------------------------------------------------------------------------#


#' Estimating the propensity score model using splines for replication Brookhart et al. (2006)
#'
#' @description This approach for using the PS to estimate exposure effects is used as one of the approaches in experiment 1 and experiment 2 (the other approach is based on subclassification). The sensitivity analyses in experiment 1 apply this approach only.
#'
#' @param PS The propensity score, as estimated by the function estimate_ps().
#' @param data The dataset containing the exposure variable A, outcome variable Y and covariates.
#'
#' @return Returns the exposure effect, estimated in an outcome model that is adjusted for the propensity score spline. Also stores warning for fitting the outcome model. Sometimes the model is overspecified; in that case the exposure effect is set to 'NA'.
#'
#' @import splines
#' @import simsalapar
estimate_effect_spline <- function(PS, data){
  # Estimate propensity score, using cubic splines with three interior knots
  PS_spline     <- splines::bs(PS, knots = quantile(PS, probs = c(0.25, 0.5, 0.75)))

  # Estimate outcome Poisson model using propensity score spline
  outcome_model <- simsalapar::tryCatch.W.E(
                    glm(Y ~ A + PS_spline, data = data, family = "poisson"))

  # Store exposure effect: store NA when warning occurred (model overspecified)
  effect_A      <- ifelse(is.null(outcome_model$value$coefficients["A"]), NA, outcome_model$value$coefficients["A"])

  if(is.null(outcome_model$warning)){outcome_model$warning <- NA}

  return(list(effect_A = effect_A, warning = outcome_model$warning))
}

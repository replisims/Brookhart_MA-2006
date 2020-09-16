#------------------------------------------------------------------------------#
# Replication of Brookhart M.A. et al. (2006)
# Variable Selection for Propensity Score Models.
# American Journal of Epidemiology, 163(12), 1149–1156.
#
# Replicator: K Luijken
# Co-pilot: B B L Penning de Vries
#
# Helper function to estimate outcome model including the ps based on subclassification
#------------------------------------------------------------------------------#


#' Estimating the propensity score model based on subclassification for replication Brookhart et al. (2006)
#'
#' @param PS The propensity score, as estimated by the function estimate_ps().
#' @param data The dataset containing the exposure variable A, outcome variable Y and relevant covariates.
#'
#' @return Returns the exposure effect, estimated within strata defined by quintiles of the propensity score and then averaged across strata. Note that NAs are removed in computing the mean;
#' NAs occur in subsets of the data that include non-exposed or exposed only, i.e., in separated datasets. Stores potential warnings for all models fitted in data subsets.
#'
#' @import data.table
#' @import simsalapar
estimate_effect_quintiles <- function(PS, data){
  df       <- data.frame(Y = data$Y, A = data$A, PS = PS)

  # Split propensity score into quintiles
  PS_quint <- quantile(df$PS, probs = seq(from = 0,to = 1, by = 0.2))
  # Subset dataset accordingly
  df$pscoreq <- cut(df$PS, breaks = PS_quint, labels = 1:5, include.lowest = T)

  # Create output matrices to store exposure effect in each quintile and potential warnings
  quintile_effect <- matrix(NA, nrow=5, ncol = 1)
  warnings <- data.table::data.table(matrix(NA, nrow=5, ncol = 1))
  for(i in 1:5){
    quintile <- df[df$pscoreq == i,]
    quintile_mod <- simsalapar::tryCatch.W.E(
                      glm(Y ~ A, data = quintile, family = poisson)) # geen  + quintile$PS, right? Mee eens!
    quintile_effect[i,] <- quintile_mod$value$coefficients["A"]
    warnings[i,] <- ifelse(is.null(quintile_mod$warning), NA, quintile_mod)
  }

  # Return mean exposure effect (note that NAs are removed in computing the mean!
  # NAs occur in subsets of the data that include non-exposed or exposed only,
  # i.e., in separated datasets)


  return(list(effect_A = mean(quintile_effect, na.rm = T), warning = warnings))
}

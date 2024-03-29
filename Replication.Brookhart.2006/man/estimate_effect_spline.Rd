% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_effect_spline.R
\name{estimate_effect_spline}
\alias{estimate_effect_spline}
\title{Estimating the propensity score model using splines for replication Brookhart et al. (2006)}
\usage{
estimate_effect_spline(PS, data)
}
\arguments{
\item{PS}{The propensity score, as estimated by the function estimate_ps().}

\item{data}{The dataset containing the exposure variable A, outcome variable Y and covariates.}
}
\value{
Returns the exposure effect, estimated in an outcome model that is adjusted for the propensity score spline. Also stores warning for fitting the outcome model. Sometimes the model is overspecified; in that case the exposure effect is set to 'NA'.
}
\description{
This approach for using the PS to estimate exposure effects is used as one of the approaches in experiment 1 and experiment 2 (the other approach is based on subclassification). The sensitivity analyses in experiment 1 apply this approach only.
}

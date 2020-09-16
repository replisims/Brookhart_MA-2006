#------------------------------------------------------------------------------#
# Replication of Brookhart M.A. et al. (2006)
# Variable Selection for Propensity Score Models.
# American Journal of Epidemiology, 163(12), 1149–1156.
#
# Replicator: K Luijken
# Co-pilot: B B L Penning de Vries
#
# Script to generate visualisations of results
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Helper functions ----
#------------------------------------------------------------------------------#

variance <- function(x) (1/(rep-1)) * sum((x - mean(x))^2)
mse <- function(x) mean((x - alpha4)^2)


#------------------------------------------------------------------------------#
# Figure 2, variance
#------------------------------------------------------------------------------#

create_figure2 <- function(){
  output <- main_scenarios_2
  output$alpha1 <- exp(output$alpha1)
  output$beta1 <- exp(output$beta1)
  output$variance_gamma_0 <- NA
  output$variance_gamma_1 <- NA

  for(i in 1:nrow(main_scenarios_2)){
    results <- readRDS(paste0("./analysis/data/raw_data/experiment_2_main/S",i,".rds"))
    results <- data.frame(results)

    output[i,"variance_gamma_0"] <- variance(results[,"effect_crude"])
    output[i,"variance_gamma_1"] <- variance(results[,"effect_splines"])
  }


  input_plot <- matrix(NA, nrow = length(beta1), ncol = 5)
  colnames(input_plot) <- c("beta1","variance_gamma_0_n500","variance_gamma_0_n2500","variance_gamma_1_n500","variance_gamma_1_n2500")

  input_plot[,"beta1"] <- beta1
  input_plot[,"variance_gamma_0_n500"] <- output[output[,"nobs"]==500  & main_scenarios_2[,"alpha1"]==0,"variance_gamma_0"]
  input_plot[,"variance_gamma_1_n500"] <- output[output[,"nobs"]==500  & main_scenarios_2[,"alpha1"]==0,"variance_gamma_1"]
  input_plot[,"variance_gamma_0_n2500"] <- output[output[,"nobs"]==2500 & main_scenarios_2[,"alpha1"]==0,"variance_gamma_0"]
  input_plot[,"variance_gamma_1_n2500"] <- output[output[,"nobs"]==2500 & main_scenarios_2[,"alpha1"]==0,"variance_gamma_1"]

  png("./analysis/figures/Figure2.png",width = 480, height = 480)
  plot(x = input_plot[,"beta1"],
       y = input_plot[,"variance_gamma_0_n500"],
       type = "l",
       lty = 4,
       xlim=c(min(beta1),max(beta1)),
       ylim = c(0,0.007),
       xlab = "Association between X1 and A (relative risk)",
       ylab = "Variance of estimator")
  lines(x = input_plot[,"beta1"], y = input_plot[,"variance_gamma_0_n2500"], lty = 3)
  lines(x = input_plot[,"beta1"], y = input_plot[,"variance_gamma_1_n500"], lty = 1)
  lines(x = input_plot[,"beta1"], y = input_plot[,"variance_gamma_1_n2500"], lty = 2)
  dev.off()

}


#------------------------------------------------------------------------------#
# Figure 3, contour plot splines
#------------------------------------------------------------------------------#

create_figure3 <- function(){
  output <- main_scenarios_2
  output$alpha1 <- exp(output$alpha1)
  output$beta1 <- exp(output$beta1)
  output$mse_gamma_0 <- NA
  output$mse_gamma_1 <- NA

  for(i in 1:nrow(main_scenarios_2)){
    results <- readRDS(paste0("./analysis/data/raw_data/experiment_2_main/S",i,".rds"))
    results <- data.frame(results)

    output[i,"mse_gamma_0"] <- mse(results[,"effect_crude"])
    output[i,"mse_gamma_1"] <- mse(results[,"effect_splines"])
  }

  ### n = 500

  output_n500 <- output[output[,"nobs"] == 500,]
  mse_ratio   <- output_n500[,"mse_gamma_1"]/output_n500[,"mse_gamma_0"]
  output_n500 <- cbind(output_n500[,c("alpha1","beta1")],mse_ratio)

  mse_loess    <- loess(mse_ratio ~ alpha1 * beta1, data = output_n500)
  grid_contour <- expand.grid(alpha1 = alpha1, beta1= beta1)
  input_contour<- predict(mse_loess, newdata = grid_contour)

  plot(x = 0, y = 0,type = "n",
       xlim = c(min(alpha1), max(alpha1)),
       ylim = c(min(beta1),max(beta1)), xlab = "", ylab = "")
  contour(x = alpha1, y = beta1, z = input_contour, levels = c(0.5,0.9,1.1,1.5))

}



#------------------------------------------------------------------------------#
# Figure 4, contour plot quintiles
#------------------------------------------------------------------------------#

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
mse <- function(x, alpha4) mean((x - alpha4)^2)


#------------------------------------------------------------------------------#
# Figure 2, variance
#------------------------------------------------------------------------------#

create_figure2 <- function(output){
  output$alpha1 <- pnorm(output$alpha1*qnorm(.75))/pnorm(output$alpha1*qnorm(.25))
  output$beta1 <- pnorm(output$beta1*qnorm(.75))/pnorm(output$beta1*qnorm(.25))
  output$variance_gamma_0 <- NA
  output$variance_gamma_1 <- NA

  for(i in 1:nrow(output)){
    results <- readRDS(paste0("./analysis/data/raw_data/experiment_2_main/S",i,".rds"))
    results <- data.frame(results)

    output[i,"variance_gamma_0"] <- var(results[,"effect_crude"])
    output[i,"variance_gamma_1"] <- var(results[,"effect_splines"])
  }


  input_plot <- matrix(NA, nrow = length(unique(output$beta1)), ncol = 5)
  colnames(input_plot) <- c("beta1","variance_gamma_0_n500","variance_gamma_0_n2500","variance_gamma_1_n500","variance_gamma_1_n2500")

  input_plot[,"beta1"] <- pnorm(beta1*qnorm(.75))/pnorm(beta1*qnorm(.25))

  for(j in 1:length(unique(output$beta1))){
    input_plot[j,"variance_gamma_0_n500"] <- mean(output[output[,"nobs"]==500 & output[,"beta1"] == input_plot[j,"beta1"],"variance_gamma_0"])
    input_plot[j,"variance_gamma_1_n500"] <- mean(output[output[,"nobs"]==500 & output[,"beta1"] == input_plot[j,"beta1"],"variance_gamma_1"])
    input_plot[j,"variance_gamma_0_n2500"] <- mean(output[output[,"nobs"]==2500 & output[,"beta1"] == input_plot[j,"beta1"],"variance_gamma_0"])
    input_plot[j,"variance_gamma_1_n2500"] <- mean(output[output[,"nobs"]==2500 & output[,"beta1"] == input_plot[j,"beta1"],"variance_gamma_1"])

  }

  png("./analysis/figures/Figure2.png",width = 480, height = 480)
  plot(x = input_plot[,"beta1"],
       y = input_plot[,"variance_gamma_0_n500"],
       type = "l",
       lty = 4,
       xlim = c(min(unique(output$beta1)), 3.1),
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

create_figure3 <- function(output, alpha4){
  output$alpha1 <- pnorm(output$alpha1*qnorm(.75))/pnorm(output$alpha1*qnorm(.25))
  output$beta1 <- pnorm(output$beta1*qnorm(.75))/pnorm(output$beta1*qnorm(.25))
  output$mse_gamma_0 <- NA
  output$mse_gamma_1 <- NA

  for(i in 1:nrow(output)){
    results <- readRDS(paste0("./analysis/data/raw_data/experiment_2_main/S",i,".rds"))
    results <- data.frame(results)

    output[i,"mse_gamma_0"] <- mse(results[,"effect_crude"], alpha4 = alpha4)
    output[i,"mse_gamma_1"] <- mse(results[,"effect_splines"], alpha4 = alpha4)
  }

  ### n = 500

  output_n500    <- output[output[,"nobs"] == 500,]
  mse_ratio_n500 <- output_n500[,"mse_gamma_1"]/output_n500[,"mse_gamma_0"]
  output_n500    <- cbind(output_n500[,c("alpha1","beta1")], mse_ratio_n500)

  mse_loess_n500    <- loess(mse_ratio_n500 ~ alpha1 * beta1, data = output_n500)
  grid_contour_n500 <- expand.grid(alpha1 = pnorm(alpha1*qnorm(.75))/pnorm(alpha1*qnorm(.25)),
                                   beta1 = pnorm(beta1*qnorm(.75))/pnorm(beta1*qnorm(.25)))
  input_contour_n500<- predict(mse_loess_n500, newdata = grid_contour_n500)

  ### n = 2,500

  output_n2500    <- output[output[,"nobs"] == 2500,]
  mse_ratio_n2500 <- output_n2500[,"mse_gamma_1"]/output_n2500[,"mse_gamma_0"]
  output_n2500    <- cbind(output_n2500[,c("alpha1","beta1")], mse_ratio_n2500)

  mse_loess_n2500    <- loess(mse_ratio_n2500 ~ alpha1 * beta1, data = output_n2500)
  grid_contour_n2500 <- expand.grid(alpha1 = pnorm(alpha1*qnorm(.75))/pnorm(alpha1*qnorm(.25)),
                                    beta1 = pnorm(beta1*qnorm(.75))/pnorm(beta1*qnorm(.25)))
  input_contour_n2500<- predict(mse_loess_n2500, newdata = grid_contour_n2500)

  png("./analysis/figures/Figure3.png",width = 480, height = 480)
  par(oma = c(3,3,3,3))
  par(mar = c(0,0,0,0))
  par(mfrow = c(1,2))
  contour(x = unique(output$alpha1), y = unique(output$beta1),
          z = input_contour_n500, levels = c(0.5,0.9,1.1,1.5),axes = FALSE, ylim = c(1,3.5))
  abline(h = 3.5)
  axis(1)
  axis(2)
  axis(3, labels = F)
  axis(4, labels = F, tick = F)
  mtext("n = 500", side=3, line=-0.7, adj=0.5, cex=0.7)
  contour(x = unique(output$alpha1), y = unique(output$beta1),
          z = input_contour_n2500, levels = c(0.5,0.9,1.1,1.5), axes = FALSE, ylim = c(1,3.5))
  abline(v = c(1,3.5))
  abline(h = 3.5)
  axis(1, labels = F)
  axis(3)
  axis(4, labels = F)
  mtext("n = 2,500", side=3, line=-0.7, adj=0.5, cex=0.7)
  box("inner")
  dev.off()
}



#------------------------------------------------------------------------------#
# Figure 4, contour plot quintiles
#------------------------------------------------------------------------------#

create_figure4 <- function(output, alpha4){
  output$alpha1 <- pnorm(output$alpha1*qnorm(.75))/pnorm(output$alpha1*qnorm(.25))
  output$beta1 <- pnorm(output$beta1*qnorm(.75))/pnorm(output$beta1*qnorm(.25))
  output$mse_gamma_0 <- NA
  output$mse_gamma_1 <- NA

  for(i in 1:nrow(scenarios_2_main)){
    results <- readRDS(paste0("./analysis/data/raw_data/experiment_2_main/S",i,".rds"))
    results <- data.frame(results)

    output[i,"mse_gamma_0"] <- mse(results[,"effect_crude"], alpha4 = alpha4)
    output[i,"mse_gamma_1"] <- mse(results[,"effect_quintiles"], alpha4 = alpha4)
  }

  ### n = 500

  output_n500    <- output[output[,"nobs"] == 500,]
  mse_ratio_n500 <- output_n500[,"mse_gamma_1"]/output_n500[,"mse_gamma_0"]
  output_n500    <- cbind(output_n500[,c("alpha1","beta1")],mse_ratio_n500)

  mse_loess_n500    <- loess(mse_ratio_n500 ~ alpha1 * beta1, data = output_n500)
  grid_contour_n500 <- expand.grid(alpha1 = pnorm(alpha1*qnorm(.75))/pnorm(alpha1*qnorm(.25)),
                                   beta1 = pnorm(beta1*qnorm(.75))/pnorm(beta1*qnorm(.25)))
  input_contour_n500<- predict(mse_loess_n500, newdata = grid_contour_n500)

  ### n = 2,500

  output_n2500    <- output[output[,"nobs"] == 2500,]
  mse_ratio_n2500 <- output_n2500[,"mse_gamma_1"]/output_n2500[,"mse_gamma_0"]
  output_n2500    <- cbind(output_n2500[,c("alpha1","beta1")],mse_ratio_n2500)

  mse_loess_n2500    <- loess(mse_ratio_n2500 ~ alpha1 * beta1, data = output_n2500)
  grid_contour_n2500 <- expand.grid(alpha1 = pnorm(alpha1*qnorm(.75))/pnorm(alpha1*qnorm(.25)),
                                   beta1 = pnorm(beta1*qnorm(.75))/pnorm(beta1*qnorm(.25)))
  input_contour_n2500<- predict(mse_loess_n2500, newdata = grid_contour_n2500)

  png("./analysis/figures/Figure4.png",width = 480, height = 480)
  par(oma = c(3,3,3,3))
  par(mar = c(0,0,0,0))
  par(mfrow = c(1,2))
  contour(x = unique(output$alpha1), y = unique(output$beta1), ylim = c(1,3.5),
          z = input_contour_n500, levels = c(0.5,0.9,1.1,1.5),axes = FALSE)
  abline(h = 3.5)
  axis(1)
  axis(2)
  axis(3, labels = F)
  axis(4, labels = F, tick = F)
  mtext("n = 500", side=3, line=-0.7, adj=0.5, cex=0.7)
  contour(x = unique(output$alpha1), y = unique(output$beta1), ylim = c(1,3.5),
          z = input_contour_n2500, levels = c(0.5,0.9,1.1,1.5), axes = FALSE)
  abline(v = c(1,3.5))
  abline(h = 3.5)
  axis(1, labels = F)
  axis(3)
  axis(4, labels = F)
  mtext("n = 2,500", side=3, line=-0.7, adj=0.5, cex=0.7)
  box("inner")
  dev.off()
}



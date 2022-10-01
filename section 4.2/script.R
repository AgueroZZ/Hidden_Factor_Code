library(SPCompute)
library(tidyverse)
library(progress)
library(latex2exp)


convert_preva_to_intercept <- function(parameters, mode, covariate, seed = 123, B = 10000){
  preva <- parameters$preva
  betaG <- parameters$betaG
  if(covariate == "binary"){

    ComputeEgivenG <- function(gamma0,gammaG,G, E = 1){
      PEG <- (exp(gamma0 + gammaG * G)^E)/(1+exp(gamma0 + gammaG * G))
      PEG
    }
    ComputeYgivenGE <- function(beta0,betaG, betaE, G, E, Y = 1){
      PYGE <- (exp(beta0 + betaG * G + betaE * E)^Y)/(1 + exp(beta0 + betaG * G + betaE * E))
      PYGE
    }

    pG <- parameters$pG
    pE <- parameters$pE
    betaE <- parameters$betaE
    gammaG <- parameters$gammaG

    if(mode == "recessive"){
      pG <- pG^2
      qG <- 1- pG
      ### Solve for gamma0
      solveForgamma0 <- function(pE,gammaG, pG){
        qG <- 1 - pG
        ComputePE <- function(gamma0){
          PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) +
            ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
          PE - pE
        }
        uniroot(ComputePE, c(-8, 8), tol = (.Machine$double.eps^0.5) )$root
      }
      ### Solve for beta0
      solveForbeta0 <- function(preva, betaG, betaE, pG, pE, gammaG){
        qG <- 1 - pG
        gamma0 <- solveForgamma0(pE,gammaG, pG)
        ComputeP <- function(beta0){
          P <- ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 0)) * (qG) +
            ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 1) * (pG) + ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 1)) * (pG)
          P - preva
        }
        uniroot(ComputeP, c(-8, 8), tol = (.Machine$double.eps^0.5) )$root
      }
      gamma0 <- solveForgamma0(pE,gammaG, pG)
      beta0 <- solveForbeta0(preva, betaG, betaE, pG, pE, gammaG)
      return(list(gamma0 = gamma0, beta0 = beta0))
    }
    else if(mode == "dominant"){
      qG <- (1-pG)^2
      pG <- 1 - qG
      ### Solve for gamma0
      solveForgamma0 <- function(pE,gammaG, pG){
        qG <- 1 - pG
        ComputePE <- function(gamma0){
          PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) +
            ComputeEgivenG(gamma0,gammaG,G = 1) * (pG)
          PE - pE
        }
        uniroot(ComputePE, c(-8, 8), tol = (.Machine$double.eps^0.5) )$root
      }
      ### Solve for beta0
      solveForbeta0 <- function(preva, betaG, betaE, pG, pE, gammaG){
        qG <- 1 - pG
        gamma0 <- solveForgamma0(pE,gammaG, pG)
        ComputeP <- function(beta0){
          P <- ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 0) * (qG) + ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 0)) * (qG) +
            ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 1) * (pG) + ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 1)) * (pG)
          P - preva
        }
        uniroot(ComputeP, c(-8, 8), tol = (.Machine$double.eps^0.5) )$root
      }
      gamma0 <- solveForgamma0(pE,gammaG, pG)
      beta0 <- solveForbeta0(preva, betaG, betaE, pG, pE, gammaG)
      return(list(gamma0 = gamma0, beta0 = beta0))
    }
    else{
      solveForgamma0 <- function(pE,gammaG, pG){
        qG <- 1 - pG
        ComputePE <- function(gamma0){
          PE <- ComputeEgivenG(gamma0,gammaG,G = 0) * (qG^2) + ComputeEgivenG(gamma0,gammaG,G = 2) * (pG^2) +
            ComputeEgivenG(gamma0,gammaG,G = 1) * (2*qG*pG)
          PE - pE
        }
        uniroot(ComputePE, c(-8, 8), tol = (.Machine$double.eps^0.5) )$root
      }
      ### Solve for beta0
      solveForbeta0 <- function(preva, betaG, betaE, pG, pE, gammaG){
        qG <- 1 - pG
        gamma0 <- solveForgamma0(pE,gammaG, pG)
        ComputeP <- function(beta0){
          P <- ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 0) * (qG^2) + ComputeYgivenGE(beta0,betaG, betaE, G = 0, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 0)) * (qG^2) +
            ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 1) * (2* pG * qG) + ComputeYgivenGE(beta0,betaG, betaE, G = 1, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0,gammaG,G = 1)) * (2* pG * qG) +
            ComputeYgivenGE(beta0,betaG, betaE, G = 2, E = 1, Y = 1) * ComputeEgivenG(gamma0,gammaG,G = 2) * (pG^2) + ComputeYgivenGE(beta0,betaG, betaE, G = 2, E = 0, Y = 1) * (1 - ComputeEgivenG(gamma0, gammaG, G = 2)) * (pG^2)
          P - preva
        }
        uniroot(ComputeP, c(-8, 8), tol = (.Machine$double.eps^0.5) )$root
      }
      gamma0 <- solveForgamma0(pE,gammaG, pG)
      beta0 <- solveForbeta0(preva, betaG, betaE, pG, pE, gammaG)
      return(list(gamma0 = gamma0, beta0 = beta0))

    }
  }
  else if(covariate == "continuous"){
    pG <- parameters$pG
    qG <- 1 - pG
    muE <- parameters$muE
    sigmaE <- parameters$sigmaE
    betaE <- parameters$betaE
    gammaG <- parameters$gammaG

    if(mode == "additive"){
      gamma0 <- muE - gammaG * (2*pG*qG + 2*pG^2)
      betaG <- parameters$betaG
      betaE <- parameters$betaE
      varG <- (2*pG*qG + 4*pG^2) - (2*pG*qG + 2*pG^2)^2
      if((sigmaE^2) <= (gammaG^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
      sigmaError <- sqrt(sigmaE^2 - (gammaG^2) * varG)
      solveForbeta0 <- function(preva, betaG, betaE, pG, gammaG){
        qG <- 1 - pG
        ComputeP <- function(beta0){
          set.seed(seed)
          G <- sample(c(0,1,2), size = B, replace = T, prob = c(qG^2, 2*pG*qG, pG^2))
          E <- gamma0 + gammaG * G + rnorm(B, sd = sigmaError)
          y <- beta0 + betaG * G + betaE * E + rlogis(B)
          P <- mean(ifelse(y > 0, 1, 0))
          P - preva
        }
        uniroot(ComputeP, c(-8, 8), tol = (.Machine$double.eps^0.5) )$root
      }
      beta0 <- solveForbeta0(preva, betaG, betaE, pG, gammaG)
      return(list(gamma0 = gamma0, sigmaError = sigmaError, beta0 = beta0))
    }
    else if(mode == "dominant"){
      pG <- 2*parameters$pG*(1-parameters$pG) + parameters$pG^2
      qG <- 1 - pG
      gammaG <- parameters$gammaG
      muE <- parameters$muE
      sigmaE <- parameters$sigmaE
      gamma0 <- muE - gammaG * (pG)
      betaG <- parameters$betaG
      betaE <- parameters$betaE
      varG <- pG*qG
      if((sigmaE^2) <= (gammaG^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
      sigmaError <- sqrt(sigmaE^2 - (gammaG^2) * varG)
      solveForbeta0 <- function(preva, betaG, betaE, pG, gammaG){
        qG <- 1 - pG
        ComputeP <- function(beta0){
          set.seed(seed)
          G <- sample(c(0,1), size = B, replace = T, prob = c(qG,pG))
          E <- gamma0 + gammaG * G + rnorm(B, sd = sigmaError)
          y <- beta0 + betaG * G + betaE * E + rlogis(B)
          P <- mean(ifelse(y > 0, 1, 0))
          P - preva
        }
        uniroot(ComputeP, c(-8, 8), tol = (.Machine$double.eps^0.5) )$root
      }
      beta0 <- solveForbeta0(preva, betaG, betaE, pG, gammaG)
      return(list(gamma0 = gamma0, sigmaError = sigmaError, beta0 = beta0))
    }
    else{
      pG <- parameters$pG^2
      qG <- 1 - pG
      gammaG <- parameters$gammaG
      muE <- parameters$muE
      sigmaE <- parameters$sigmaE
      gamma0 <- muE - gammaG * (pG)
      betaG <- parameters$betaG
      betaE <- parameters$betaE
      varG <- pG*qG
      if((sigmaE^2) <= (gammaG^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
      sigmaError <- sqrt(sigmaE^2 - (gammaG^2) * varG)

      solveForbeta0 <- function(preva, betaG, betaE, pG, gammaG){
        qG <- 1 - pG
        ComputeP <- function(beta0){
          set.seed(seed)
          G <- sample(c(0,1), size = B, replace = T, prob = c(qG,pG))
          E <- gamma0 + gammaG * G + rnorm(B, sd = sigmaError)
          y <- beta0 + betaG * G + betaE * E + rlogis(B)
          P <- mean(ifelse(y > 0, 1, 0))
          P - preva
        }
        uniroot(ComputeP, c(-8, 8), tol = (.Machine$double.eps^0.5) )$root
      }
      beta0 <- solveForbeta0(preva, betaG, betaE, pG, gammaG)
      return(list(gamma0 = gamma0, sigmaError = sigmaError, beta0 = beta0))
    }

  }
  else{
    ComputeYgivenG <- function(beta0,betaG, G, Y = 1){
      PYG <- (exp(beta0 + betaG * G)^Y)/(1 + exp(beta0 + betaG * G))
      PYG
    }
    if(mode == "dominant"){
      solveForbeta0 <- function(preva, betaG, pG){
        qG <- 1 - pG
        ComputeP <- function(beta0){
          P <- ComputeYgivenG(beta0,betaG, G = 0, Y = 1) * (qG) + ComputeYgivenG(beta0,betaG, G = 1, Y = 1) * (pG)
          P - preva
        }
        uniroot(ComputeP, c(-8, 8), tol = (.Machine$double.eps^0.5) )$root
      }
      pG <- (parameters$pG^2) + 2*((1-parameters$pG) * parameters$pG)
      beta0 <- solveForbeta0(preva, betaG, pG)
      return(beta0)
    }
    else if(mode == "additive"){
      solveForbeta0 <- function(preva, betaG, pG){
        qG <- 1 - pG
        ComputeP <- function(beta0){
          P <- ComputeYgivenG(beta0,betaG, G = 0, Y = 1) * (qG^2) + ComputeYgivenG(beta0,betaG, G = 2, Y = 1) * (pG^2) + ComputeYgivenG(beta0,betaG, G = 1, Y = 1) * (2*pG*qG)
          P - preva
        }
        uniroot(ComputeP, c(-8, 8), tol = (.Machine$double.eps^0.5) )$root
      }
      beta0 <- solveForbeta0(preva, betaG, parameters$pG)
      return(beta0)
    }
    else{
      solveForbeta0 <- function(preva, betaG, pG){
        qG <- 1 - pG
        ComputeP <- function(beta0){
          P <- ComputeYgivenG(beta0,betaG, G = 0, Y = 1) * (qG) + ComputeYgivenG(beta0,betaG, G = 1, Y = 1) * (pG)
          P - preva
        }
        uniroot(ComputeP, c(-8, 8), tol = (.Machine$double.eps^0.5) )$root
      }
      pG <- (parameters$pG^2)
      beta0 <- solveForbeta0(preva, betaG, pG)
      return(beta0)
    }

  }
}
# setwd("~/Documents/GitHub/Stats-Gene/Bioinfo_Project/Simulations/sim2/")



###### Case 1: Just G alone

#### Study the variability of estimated power in relation to B
set.seed(123)
parameters <- list(preva = 0.2, pG = 0.1, betaG = log(1.5))
others_para <- convert_preva_to_intercept(parameters, mode = "additive", covariate = "none")
beta0 <- convert_preva_to_intercept(parameters, mode = "additive", covariate = "none", seed = 2)
exp(beta0)/(1+exp(beta0))
n <- 1000

compute_SE_at_one_point <- function(B, rep_num = 100, show_result = T){
  powers <- numeric(rep_num)
  time <- numeric(rep_num)
  if(show_result){
    pb <- progress_bar$new(total = rep_num)
    for (i in 1:rep_num) {
      begin_time <- Sys.time()
      powers[i] <- SPCompute:::Compute_Power_Sim(parameters, n = n, mode = "additive", covariate = "none", response = "binary", seed = i, B = B)
      end_time <- Sys.time()
      time[i] <- as.numeric(end_time - begin_time)
      pb$tick()
    }
  }
  else{
    for (i in 1:rep_num) {
      begin_time <- Sys.time()
      powers[i] <- SPCompute:::Compute_Power_Sim(parameters, n = n, mode = "additive", covariate = "none", response = "binary", seed = i, B = B)
      end_time <- Sys.time()
      time[i] <- as.numeric(end_time - begin_time)
    }
  }
  list(SE = sd(powers), Time = mean(time))
}

set.seed(123)
B_vec <- seq(1000, 20000, by = 500)
compute_SE <- function(B_vec){
  SE_vec <- numeric(length = length(B_vec))
  Time_vec <- numeric(length = length(B_vec))
  pb <- progress_bar$new(total = length(SE_vec))
  for (i in 1:length(SE_vec)) {
    res <- compute_SE_at_one_point(B_vec[i], rep_num = 100, show_result = F)
    SE_vec[i] <- res$SE
    Time_vec[i] <- res$Time
    pb$tick()
  }
  return(data.frame(B = B_vec, SE = SE_vec, Times = Time_vec))
}


result <- compute_SE(B_vec)
pdf(file = "sim2_SEplot_case1.pdf", height = 4, width = 5)
result %>% ggplot() + geom_point(aes(x = B, y = SE)) + geom_smooth(aes(x = B, y = SE)) + theme_classic()
dev.off()
pdf(file = "sim2_logSEplot_case1.pdf", height = 4, width = 5)
result %>% mutate(logB = log(B), logSE = log(SE)) %>% ggplot() + geom_point(aes(x = logB, y = logSE), color = 2) + geom_smooth(aes(x = logB, y = logSE), method = 'lm') + theme_classic()
dev.off()
pdf(file = "sim2_Timeplot_case1.pdf", height = 4, width = 5)
result %>% ggplot() + geom_point(aes(x = B, y = Times)) + geom_smooth(aes(x = B, y = Times)) + theme_classic()
dev.off()

summary(lm(Times~B, data = result)) #1.124e-05
summary(lm(log(I(SE))~log(I(B)), data = result)) # -0.469155



#### Study the distribution of estimated power at a specific B value:
compute_Powers_vec <- function(B, rep_num = 100, show_result = T){
  powers <- numeric(rep_num)
  time <- numeric(rep_num)
  if(show_result){
    pb <- progress_bar$new(total = rep_num)
    for (i in 1:rep_num) {
      begin_time <- Sys.time()
      powers[i] <- SPCompute:::Compute_Power_Sim(parameters, n = n, mode = "additive", covariate = "none", response = "binary", seed = i, B = B)
      end_time <- Sys.time()
      time[i] <- as.numeric(end_time - begin_time)
      pb$tick()
    }
  }
  else{
    for (i in 1:rep_num) {
      begin_time <- Sys.time()
      powers[i] <- SPCompute:::Compute_Power_Sim(parameters, n = n, mode = "additive", covariate = "none", response = "binary", seed = i, B = B)
      end_time <- Sys.time()
      time[i] <- as.numeric(end_time - begin_time)
    }
  }
  data.frame(Powers = powers, Times = time)
}

pdf(file = "sim2_Power100plot_case1.pdf", height = 4, width = 5)
powers_vec_100 <- compute_Powers_vec(B = 100, rep_num = 1000)
powers_vec_100 %>% ggplot() + geom_density(aes(x = Powers), bw = 0.04) + geom_histogram(aes(x = Powers, y = ..density..), bins = 10, alpha = 0.3) + theme_classic()
dev.off()

pdf(file = "sim2_Power10000plot_case1.pdf", height = 4, width = 5)
powers_vec_10000 <- compute_Powers_vec(B = 10000, rep_num = 1000)
powers_vec_10000 %>% ggplot() + geom_density(aes(x = Powers), bw = 0.006) + geom_histogram(aes(x = Powers, y = ..density..), bins = 10, alpha = 0.3) + theme_classic()
dev.off()

pdf(file = "sim2_Power20000plot_case1.pdf", height = 4, width = 5)
powers_vec_20000 <- compute_Powers_vec(B = 20000, rep_num = 1000)
powers_vec_20000 %>% ggplot() + geom_density(aes(x = Powers), bw = 0.006) + geom_histogram(aes(x = Powers, y = ..density..), bins = 10, alpha = 0.3) + theme_classic()
dev.off()


###### Case 2: Binary E

#### Study the variability of estimated power in relation to B
set.seed(123)
parameters <- list(preva = 0.2, pE = 0.3, pG = 0.1, betaG = log(1.5), betaE = log(2.5), gammaG = log(0.2))
others_para <- convert_preva_to_intercept(parameters, mode = "additive", covariate = "binary")
beta0 <- convert_preva_to_intercept(parameters, mode = "additive", covariate = "binary")$beta0
exp(beta0)/(1+exp(beta0))
n <- 1000

compute_SE_at_one_point <- function(B, rep_num = 100, show_result = T){
  powers <- numeric(rep_num)
  time <- numeric(rep_num)
  if(show_result){
    pb <- progress_bar$new(total = rep_num)
    for (i in 1:rep_num) {
      begin_time <- Sys.time()
      powers[i] <- SPCompute:::Compute_Power_Sim(parameters, n = n, mode = "additive", covariate = "binary", response = "binary", seed = i, B = B)
      end_time <- Sys.time()
      time[i] <- as.numeric(end_time - begin_time)
      pb$tick()
    }
  }
  else{
    for (i in 1:rep_num) {
      begin_time <- Sys.time()
      powers[i] <- SPCompute:::Compute_Power_Sim(parameters, n = n, mode = "additive", covariate = "binary", response = "binary", seed = i, B = B)
      end_time <- Sys.time()
      time[i] <- as.numeric(end_time - begin_time)
    }
  }
  list(SE = sd(powers), Time = mean(time))
}

set.seed(123)
B_vec <- seq(1000, 20000, by = 500)
compute_SE <- function(B_vec){
  SE_vec <- numeric(length = length(B_vec))
  Time_vec <- numeric(length = length(B_vec))
  pb <- progress_bar$new(total = length(SE_vec))
  for (i in 1:length(SE_vec)) {
    res <- compute_SE_at_one_point(B_vec[i], rep_num = 100, show_result = F)
    SE_vec[i] <- res$SE
    Time_vec[i] <- res$Time
    pb$tick()
  }
  return(data.frame(B = B_vec, SE = SE_vec, Times = Time_vec))
}





result <- compute_SE(B_vec)
pdf(file = "sim2_SEplot_case2.pdf", height = 4, width = 5)
result %>% ggplot() + geom_point(aes(x = B, y = SE)) + geom_smooth(aes(x = B, y = SE)) + theme_classic()
dev.off()
pdf(file = "sim2_log10SEplot_case2.pdf", height = 4, width = 5)
result %>% mutate(logB = log10(B), logSE = log10(SE)) %>% ggplot() + geom_point(aes(x = logB, y = logSE), color = 2) + 
  geom_smooth(aes(x = logB, y = logSE), method = 'lm') + theme_classic() + ylab(TeX("$\\log_{10}SE$")) +  xlab(TeX("$\\log_{10}B$"))
dev.off()
pdf(file = "sim2_log10SEplot_binaryE.pdf", height = 4, width = 5)
result %>% mutate(logB = log10(B), logSE = log10(SE)) %>% ggplot() + geom_point(aes(x = logB, y = logSE), color = 2) + 
  geom_smooth(aes(x = logB, y = logSE), method = 'lm') + theme_classic() + ylab(TeX("$\\log_{10}SE$")) +  xlab(TeX("$\\log_{10}B$"))
dev.off()
pdf(file = "sim2_Timeplot_case2.pdf", height = 4, width = 5)
result %>% ggplot() + geom_point(aes(x = B, y = Times)) + geom_smooth(aes(x = B, y = Times)) + theme_classic()
dev.off()

summary(lm(Times~B, data = result)) #1.581e-05
summary(lm(log(I(SE))~log(I(B)), data = result)) # -0.488492





#### Study the distribution of estimated power at a specific B value:
compute_Powers_vec <- function(B, rep_num = 100, show_result = T){
  powers <- numeric(rep_num)
  time <- numeric(rep_num)
  if(show_result){
    pb <- progress_bar$new(total = rep_num)
    for (i in 1:rep_num) {
      begin_time <- Sys.time()
      powers[i] <- SPCompute:::Compute_Power_Sim(parameters, n = n, mode = "additive", covariate = "binary", response = "binary", seed = i, B = B)
      end_time <- Sys.time()
      time[i] <- as.numeric(end_time - begin_time)
      pb$tick()
    }
  }
  else{
    for (i in 1:rep_num) {
      begin_time <- Sys.time()
      powers[i] <- SPCompute:::Compute_Power_Sim(parameters, n = n, mode = "additive", covariate = "none", response = "binary", seed = i, B = B)
      end_time <- Sys.time()
      time[i] <- as.numeric(end_time - begin_time)
    }
  }
  data.frame(Powers = powers, Times = time)
}

pdf(file = "sim2_Power100plot_case2.pdf", height = 4, width = 5)
powers_vec_100 <- compute_Powers_vec(B = 100, rep_num = 1000)
powers_vec_100 %>% ggplot() + geom_density(aes(x = Powers), bw = 0.04) + geom_histogram(aes(x = Powers, y = ..density..), bins = 10, alpha = 0.3) + theme_classic()
dev.off()

pdf(file = "sim2_Power10000plot_case2.pdf", height = 4, width = 5)
powers_vec_10000 <- compute_Powers_vec(B = 10000, rep_num = 1000)
powers_vec_10000 %>% ggplot() + geom_density(aes(x = Powers), bw = 0.006) + geom_histogram(aes(x = Powers, y = ..density..), bins = 10, alpha = 0.3) + theme_classic()
dev.off()

pdf(file = "sim2_Power20000plot_case2.pdf", height = 4, width = 5)
powers_vec_20000 <- compute_Powers_vec(B = 20000, rep_num = 1000)
powers_vec_20000 %>% ggplot() + geom_density(aes(x = Powers), bw = 0.006) + geom_histogram(aes(x = Powers, y = ..density..), bins = 10, alpha = 0.3) + theme_classic()
dev.off()


############################ On the selection of method: ##############################
#### Fix B to be 10000, but vary the sample size.

#### Case 1:
set.seed(123)
parameters <- list(preva = 0.2, pG = 0.1, betaG = log(1.5))
others_para <- convert_preva_to_intercept(parameters, mode = "additive", covariate = "none")
beta0 <- convert_preva_to_intercept(parameters, mode = "additive", covariate = "none", seed = 2)

run_time_comparison <- function(nvec, each = 10){
  k <- length(nvec)
  time_m1 <- c()
  time_m2 <- c()
  pb <- progress_bar$new(total = k)
  for (i in 1:k) {
    time_m1_replicate <- c()
    time_m2_replicate <- c()
    for (j in 1:each) {
      begin <- Sys.time()
      p1 <- SPCompute:::Compute_Power_Sim(parameters, n = nvec[i], mode = "additive", covariate = "none", response = "binary")
      end <- Sys.time()
      time_m1_replicate[j] <- as.numeric(end - begin)
      begin <- Sys.time()
      p2 <- SPCompute:::Compute_Power_Expanded(parameters, n = nvec[i], mode = "additive", covariate = "none", response = "binary")
      end <- Sys.time()
      time_m2_replicate[j] <- as.numeric(end - begin)
    }
    time_m1[i] <- mean(time_m1_replicate)
    time_m2[i] <- mean(time_m2_replicate)
    pb$tick()
  }
  data.frame(Method1 = time_m1, Method2 = time_m2)
}

nvec <- seq(1000, 100000, by = 1000)
times_result <- run_time_comparison(nvec)
times_result$n <- nvec


pdf(file = "sim2_Select_case1.pdf", height = 4, width = 5)
times_result %>% pivot_longer(cols = Method1:Method2, names_to = "Method", values_to = "Times") %>% 
  mutate(Method = as.factor(Method)) %>% 
  ggplot() + geom_line(aes(x = n, y = Times, color = Method)) + theme_classic()
dev.off()



#### Case 2:
set.seed(123)
parameters <- list(preva = 0.2, pE = 0.3, pG = 0.1, betaG = log(1.5), betaE = log(2.5), gammaG = log(0.2))
others_para <- convert_preva_to_intercept(parameters, mode = "additive", covariate = "binary")
beta0 <- convert_preva_to_intercept(parameters, mode = "additive", covariate = "binary", seed = 2)$beta0

run_time_comparison <- function(nvec, each = 10){
  k <- length(nvec)
  time_m1 <- c()
  time_m2 <- c()
  pb <- progress_bar$new(total = k)
  for (i in 1:k) {
    time_m1_replicate <- c()
    time_m2_replicate <- c()
    for (j in 1:each) {
      begin <- Sys.time()
      p1 <- SPCompute:::Compute_Power_Sim(parameters, n = nvec[i], mode = "additive", covariate = "binary", response = "binary")
      end <- Sys.time()
      time_m1_replicate[j] <- as.numeric(end - begin)
      begin <- Sys.time()
      p2 <- SPCompute:::Compute_Power_Expanded(parameters, n = nvec[i], mode = "additive", covariate = "binary", response = "binary")
      end <- Sys.time()
      time_m2_replicate[j] <- as.numeric(end - begin)
    }
    time_m1[i] <- mean(time_m1_replicate)
    time_m2[i] <- mean(time_m2_replicate)
    pb$tick()
  }
  data.frame(Proposed1 = time_m1, Proposed2 = time_m2)
}

nvec <- seq(1000, 100000, by = 1000)
times_result <- run_time_comparison(nvec)
times_result$n <- nvec


pdf(file = "sim2_Select_binaryE_alternative.pdf", height = 4, width = 5)
times_result %>% pivot_longer(cols = Proposed1:Proposed2, names_to = "Method", values_to = "Times") %>% 
  mutate(Method = as.factor(Method)) %>% 
  ggplot() + geom_line(aes(x = n, y = Times, color = Method)) + theme_classic() + ylab("Seconds Per Computation") + 
  scale_colour_manual(labels = c("Semi-Simulation", "Representative Data"), values = c("red", "blue"))
dev.off()


pdf(file = "sim2_Select_binaryE.pdf", height = 4, width = 5)
times_result %>% pivot_longer(cols = Proposed1:Proposed2, names_to = "Method", values_to = "Times") %>% 
  mutate(Method = as.factor(Method)) %>% 
  ggplot() + geom_line(aes(x = n, y = Times, color = Method)) + theme_classic() + ylab("Seconds Per Computation")
dev.off()




































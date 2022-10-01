library(SPCompute)
library(tidyverse)
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
compute_oracle_size <- function(Target_Power = 0.8, alpha = 0.05, OR, setting){
  if(setting == 1){
    ### Generate prospective data
    n <- seq(from = 500, to = 41000, by = 250)
    parameters <- list(preva = 0.2, pE = 0.3, pG = 0.1, betaG = log(OR), betaE = log(2.5), gammaG = log(0.2))
    others_para <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "binary")
    emp_powers <- c()
    for (i in 1:length(n)) {
      correct <- c()
      for (j in 1:1000) {
        k <- n[i]
        G <- sample(c(0,1), size = k, replace = T, prob = c((1-parameters$pG)^2, (1 - ((1-parameters$pG)^2))))
        E <- others_para$gamma0 + parameters$gammaG * G + rlogis(k)
        E <- ifelse(E>=0, 1, 0)
        y <- others_para$beta0 + parameters$betaG*G + parameters$betaE*E + rlogis(k)
        y <- ifelse(y>=0, 1, 0)
        mod <- glm(y ~ G + E, family = binomial(link = "logit"))
        correct[j] <- summary(mod)$coefficients[2,4] <= alpha
      }
      emp_powers[i] <- mean(correct)
      if(emp_powers[i] >= Target_Power) return(k)
    }
  }
  else if(setting == 2){
    ### Generate prospective data
    parameters <- list(preva = 0.2, pG = 0.1, muE = 0, sigmaE = 1, betaG = log(OR), betaE = log(2.5), gammaG = log(0.5))
    others_para <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "continuous")
    beta0 <- others_para$beta0
    n <- seq(from = 500, to = 55000, by = 250)
    emp_powers <- c()
    for (i in 1:length(n)) {
      correct <- c()
      for (j in 1:1000) {
        k <- n[i]
        G <- sample(c(0,1), size = k, replace = T, prob = c((1-parameters$pG)^2, (1 - ((1-parameters$pG)^2))))
        E <- others_para$gamma0 + parameters$gammaG * G + rnorm(n = k, sd = others_para$sigmaError)
        y <- others_para$beta0 + parameters$betaG*G + parameters$betaE*E + rlogis(k)
        y <- ifelse(y>0, 1, 0)
        mod <- glm(y~G+E, family = binomial(link = "logit"))
        correct[j] <- summary(mod)$coefficients[2,4] <= alpha
      }
      emp_powers[i] <- mean(correct)
      if(emp_powers[i] >= Target_Power) return(k)
    }
  }
  else{
    parameters <- list(preva = 0.5, pG = 0.1, betaG = log(OR))
    others_para <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "none")
    beta0 <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "none", seed = 2)
    N <- c(seq(from = 100, to = 1100, by = 40), seq(from = 1300, to = 1600, by = 20), seq(from = 3000, to = 3200, by = 20), seq(from = 10000, to = 12000, by = 50))
    n <- 2*N
    ### Generate retrospective data, compute the power curve
    emp_powers <- c()
    for (i in 1:length(n)) {
      correct <- c()
      for (j in 1:1000) {
        if(n[i] >= 2000){
          B <- 80000
        }
        else{
          B <- 10000
        }
        G <- sample(c(0,1), size = B, replace = T, prob = c((1-parameters$pG)^2, (1 - ((1-parameters$pG)^2))))
        y <- beta0 + parameters$betaG*G + rlogis(B)
        y <- ifelse(y>=0, 1, 0)
        data <- data.frame(G, y)
        case_data <- data %>% filter(y == 1) %>% sample_n(N[i], replace = F)
        control_data <- data %>% filter(y == 0) %>% sample_n((n[i] - N[i]), replace = F)
        sampled_data <- rbind(case_data, control_data)
        mod <- glm(y~G, family = binomial(link = "logit"), data = sampled_data)
        correct[j] <- summary(mod)$coefficients[2,4] <= alpha
      }
      emp_powers[i] <- mean(correct)
      if(emp_powers[i] >= Target_Power) return(n[i])
    }
  }
  print(emp_powers)
  return(message("No appropriate sample size found!"))
}
compute_oracle_size_2 <- function(Target_Power = 0.8, alpha = 0.05, OR, setting){
  if(setting == 1){
    ### Generate prospective data
    nstart <- 100
    parameters <- list(preva = 0.2, pE = 0.3, pG = 0.1, betaG = log(OR), betaE = log(2.5), gammaG = log(0.2))
    others_para <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "binary")
    emp_powers <- c()
    n <- c()
    find_ans <- FALSE
    while(find_ans == FALSE){
      correct <- c()
      for (j in 1:1000) {
        step_size <- round(0.1*nstart)
        k <- nstart + step_size
        G <- sample(c(0,1), size = k, replace = T, prob = c((1-parameters$pG)^2, (1 - ((1-parameters$pG)^2))))
        E <- others_para$gamma0 + parameters$gammaG * G + rlogis(k)
        E <- ifelse(E>=0, 1, 0)
        y <- others_para$beta0 + parameters$betaG*G + parameters$betaE*E + rlogis(k)
        y <- ifelse(y>=0, 1, 0)
        mod <- glm(y ~ G + E, family = binomial(link = "logit"))
        correct[j] <- summary(mod)$coefficients[2,4] <= alpha
      }
      emp_powers <- c(emp_powers, mean(correct))
      indx <- length(emp_powers)
      n <- c(n, k)
      if(emp_powers[indx] >= (Target_Power-0.05) & emp_powers[indx] <= (Target_Power + 0.05)) {find_ans <- TRUE}
      else if(emp_powers[indx] < (Target_Power-0.3)){
        nstart <- 4*nstart
      }
      else if(emp_powers[indx] > (Target_Power+0.3)){
        nstart <- round(nstart/4)
      }
    }
  }
  else if(setting == 2){
    ### Generate prospective data
    parameters <- list(preva = 0.2, pG = 0.1, muE = 0, sigmaE = 1, betaG = log(OR), betaE = log(2.5), gammaG = log(0.5))
    others_para <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "continuous")
    beta0 <- others_para$beta0
    nstart <- 100
    emp_powers <- c()
    n <- c()
    find_ans <- FALSE
    while(find_ans == FALSE){
      correct <- c()
      for (j in 1:1000) {
        step_size <- round(0.1*nstart)
        k <- nstart + step_size
        G <- sample(c(0,1), size = k, replace = T, prob = c((1-parameters$pG)^2, (1 - ((1-parameters$pG)^2))))
        E <- others_para$gamma0 + parameters$gammaG * G + rnorm(n = k, sd = others_para$sigmaError)
        y <- others_para$beta0 + parameters$betaG*G + parameters$betaE*E + rlogis(k)
        y <- ifelse(y>0, 1, 0)
        mod <- glm(y~G+E, family = binomial(link = "logit"))
        correct[j] <- summary(mod)$coefficients[2,4] <= alpha
      }
      emp_powers <- c(emp_powers, mean(correct))
      indx <- length(emp_powers)
      n <- c(n, k)
      if(emp_powers[indx] >= (Target_Power-0.05) & emp_powers[indx] <= (Target_Power + 0.05)) {find_ans <- TRUE}
      else if(emp_powers[indx] < (Target_Power-0.3)){
        nstart <- 4*nstart
      }
      else if(emp_powers[indx] > (Target_Power+0.3)){
        nstart <- round(nstart/4)
      }
    }
  }
  else{
    parameters <- list(preva = 0.2, pG = 0.1, betaG = log(OR))
    Nstart <- 50
    ### Generate retrospective data, compute the power curve
    emp_powers <- c()
    n <- c()
    find_ans <- FALSE
    while(find_ans == FALSE){
      correct <- c()
      for (j in 1:1000) {
        step_size <- round(0.1*Nstart)
        N <- Nstart + step_size
        k <- 5*N
        if(k >= 5000){
          B <- 60000
        }
        else{
          B <- 10000
        }
        G <- sample(c(0,1), size = B, replace = T, prob = c((1-parameters$pG)^2, (1 - ((1-parameters$pG)^2))))
        y <- -2 + parameters$betaG*G + rlogis(B)
        y <- ifelse(y>=0, 1, 0)
        data <- data.frame(G, y)
        case_data <- data %>% filter(y == 1) %>% sample_n(N, replace = F)
        control_data <- data %>% filter(y == 0) %>% sample_n((k - N), replace = F)
        sampled_data <- rbind(case_data, control_data)
        mod <- glm(y~G, family = binomial(link = "logit"), data = sampled_data)
        correct[j] <- summary(mod)$coefficients[2,4] <= alpha
      }
      emp_powers <- c(emp_powers, mean(correct))
      indx <- length(emp_powers)
      n <- c(n, k)
      if(emp_powers[indx] >= (Target_Power-0.05) & emp_powers[indx] <= (Target_Power + 0.05)) {find_ans <- TRUE}
      else if(emp_powers[indx] < (Target_Power-0.3)){
        Nstart <- 4*Nstart
      }
      else if(emp_powers[indx] > (Target_Power+0.3)){
        Nstart <- round(Nstart/4)
      }
    }
  }
  return(k)
}
setwd("~/Documents/GitHub/Stats-Gene/Bioinfo_Project/Simulations/sim1")


########################## Comparison number 1 (Prospective Sampling): ##################################
### We study two cases:
### 1. When sample sizes are given, compare the computed powers with empirical powers
### 2. When powers are given, compute and compare the computed required sample sizes
set.seed(1234)
parameters <- list(preva = 0.2, pE = 0.3, pG = 0.1, betaG = log(1.5), betaE = log(2.5), gammaG = log(0.2))
others_para <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "binary")
beta0 <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "binary")$beta0
exp(beta0)/(1+exp(beta0))
### First using the two proposed methods:
n <- seq(from = 1000, to = 3800, by = 200)
Power_Proposed1 <- c()
Power_Proposed2 <- c()
for (i in 1:length(n)) {
  Power_Proposed1[i] <- SPCompute:::Compute_Power_Sim(parameters, n = n[i], mode = "dominant", covariate = "binary", response = "binary")
  Power_Proposed2[i] <- SPCompute:::Compute_Power_Expanded(parameters, n = n[i], mode = "dominant", covariate = "binary", response = "binary")
}
### Compare the three power curves:
sim1_data_case1_givenSize <- data.frame(n = n)
sim1_data_case1_givenSize$n <- n
sim1_data_case1_givenSize$Demidenko <- c(0.5, 0.574, 0.64, 0.698, 0.748, 0.791, 0.828, 0.859, 0.885, 0.906, 0.924, 0.939, 0.951, 0.961, 0.969)
sim1_data_case1_givenSize$Quanto <- c(0.557, 0.6349, 0.7018, 0.7585, 0.8059, 0.8451, 0.8771, 0.9031, 0.9240, 0.9407, 0.9539, 0.9644, 0.9725, 0.9789, 0.9839)
sim1_data_case1_givenSize$Proposed1 <- Power_Proposed1
sim1_data_case1_givenSize$Proposed2 <- Power_Proposed2
### Compute empirical powers: (each from 1000 replications)
emp_powers <- c()
for (i in 1:length(n)) {
  correct <- c()
  for (j in 1:1000) {
    k <- n[i]
    G <- sample(c(0,1), size = k, replace = T, prob = c((1-parameters$pG)^2, (1 - ((1-parameters$pG)^2))))
    E <- others_para$gamma0 + parameters$gammaG * G + rlogis(k)
    E <- ifelse(E>=0, 1, 0)
    y <- others_para$beta0 + parameters$betaG*G + parameters$betaE*E + rlogis(k)
    y <- ifelse(y>=0, 1, 0)
    mod <- glm(y ~ G + E, family = binomial(link = "logit"))
    correct[j] <- summary(mod)$coefficients[2,4] <= 0.05
  }
  emp_powers[i] <- mean(correct)
}
sim1_data_case1_givenSize$Empirical <- emp_powers
sim1_data_case1_givenSize <- round(sim1_data_case1_givenSize, 3)
save(file = "sim1_data_case1_givenSize.rda", sim1_data_case1_givenSize)
pdf(file = "sim1_case1_power", height = 4, width = 5)
plot(emp_powers~n, type = 'l', ylab = "Estimated Power", xlab = "Sample Size", lty = "solid")
lines(Power_Proposed1~n, col = "red", lty = "dashed")
lines(Power_Proposed2~n, col = "blue", lty = "dotted")
lines(sim1_data_case1_givenSize$Demidenko~n, col = "green", lty = "dotdash")
lines(sim1_data_case1_givenSize$Quanto~n, col = "purple", lty = "longdash")
legend(x = 2700, y = 0.75,   # Coordinates (x also accepts keywords)
       legend = c("Oracle", "P1.SS", "P2.RD", "Demidenko", "Quanto"), # Vector with the name of each group
       col = c("black", "red", "blue", "green", "purple"), bty = "n", cex = 0.8,
       lty = c("solid", "dashed", "dotted", "dotdash", "longdash")
)
dev.off()
methods <- c("Pro1", "Pro2", "Demidenko", "Quanto")
# CR <- c(mean(upper >= Power_Proposed1 & Power_Proposed1 >= lower), mean(upper >= Power_Proposed2 & Power_Proposed2 >= lower), mean(upper >= sim1_data_case1_givenSize$Demidenko & sim1_data_case1_givenSize$Demidenko >= lower), mean(upper >= sim1_data_case1_givenSize$Quanto & sim1_data_case1_givenSize$Quanto >= lower))
MAE <- c(mean(abs(emp_powers - Power_Proposed1)), mean(abs(emp_powers - Power_Proposed2)), mean(abs(emp_powers -sim1_data_case1_givenSize$Demidenko)), mean(abs(emp_powers -sim1_data_case1_givenSize$Quanto)))
MaxAE <- c(max(abs(emp_powers - Power_Proposed1)), max(abs(emp_powers - Power_Proposed2)), max(abs(emp_powers -sim1_data_case1_givenSize$Demidenko)), max(abs(emp_powers -sim1_data_case1_givenSize$Quanto)))
summary_data <- data.frame(methods = methods, MAE = MAE, MaxAE = MaxAE)
summary_data
### Compare the three sample size curves (fix the power to be at 0.8, vary the effect size or MAF):
sim1_data_case1_givenPower <- data.frame(RG = seq(1.1,2.5, by = 0.1))
n_Proposed1 <- c()
n_Proposed2 <- c()
set.seed(123)
Pry <- c()
for (i in 1:length(sim1_data_case1_givenPower$RG)) {
  parameters$betaG <- log(sim1_data_case1_givenPower$RG[i])
  others_para <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "binary")
  beta0 <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "binary", seed = 2)$beta0
  Pry[i] <- exp(beta0)/(1+exp(beta0))
  n_Proposed1[i] <- SPCompute:::Compute_Size_Sim(parameters, PowerAim = 0.8, mode = "dominant", covariate = "binary", response = "binary")
}
for (i in 1:length(sim1_data_case1_givenPower$RG)) {
  parameters$betaG <- log(sim1_data_case1_givenPower$RG[i])
  others_para <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "binary")
  beta0 <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "binary", seed = 2)$beta0
  Pry[i] <- exp(beta0)/(1+exp(beta0))
  n_Proposed2[i] <- SPCompute:::Compute_Size_Expanded(parameters, PowerAim = 0.8, mode = "dominant", covariate = "binary", response = "binary", lower.lim.n = 300)
}
sim1_data_case1_givenPower$Pro1 <- ceiling(n_Proposed1)
sim1_data_case1_givenPower$Pro2 <- ceiling(n_Proposed2)
sim1_data_case1_givenPower$Quanto <- 5 * c(6858, 1837, 872, 522, 355, 261, 202, 163, 136, 115, 100, 88, 78, 71, 64)
sim1_data_case1_givenPower$Demidenko <- c(40899, 10742, 5012, 3098, 2075, 1507, 1157, 924, 762, 643, 578, 505, 447, 400, 361)
start_time <- Sys.time()
sim1_data_case1_givenPower$oracle <- c()
for (i in 1:length(sim1_data_case1_givenPower$RG)) {
  set.seed(123)
  sim1_data_case1_givenPower$oracle[i] <- compute_oracle_size(OR = sim1_data_case1_givenPower$RG[i], setting = 1)
}
end_time <- Sys.time()
end_time - start_time
save(file = "sim1_data_case1_givenPower.rda", sim1_data_case1_givenPower)
pdf(file = "sim1_case1_ss", height = 4, width = 5)
plot(sim1_data_case1_givenPower$Pro1 ~ sim1_data_case1_givenPower$RG, type = 'l', ylab = "Estimated Sample Size", col = "red", xlab = "OR", lty = "dashed")
lines(sim1_data_case1_givenPower$oracle ~ sim1_data_case1_givenPower$RG, col = "black", lty = "solid")
lines(sim1_data_case1_givenPower$Pro2 ~ sim1_data_case1_givenPower$RG, col = "blue", lty = "dotted")
lines(sim1_data_case1_givenPower$Demidenko ~ sim1_data_case1_givenPower$RG, col = "green", lty = "dotdash")
lines(sim1_data_case1_givenPower$Quanto ~ sim1_data_case1_givenPower$RG, col = "purple", lty = "longdash")
legend(x = 2, y = 25000,   # Coordinates (x also accepts keywords)
       legend = c("Oracle","P1.SS", "P2.RD", "Demidenko", "Quanto"), # Vector with the name of each group
       lty = c("solid", "dashed", "dotted", "dotdash", "longdash"),
       col = c("black", "red", "blue" , "green", "purple"), bty = "n", cex = 0.8
)
dev.off()





########################## Comparison number 2: (dominant with continuous E) #############
###############################################################################################
set.seed(123)
parameters <- list(preva = 0.2, pG = 0.1, muE = 0, sigmaE = 1, betaG = log(1.3), betaE = log(2.5), gammaG = log(0.5))
others_para <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "continuous")
beta0 <- others_para$beta0
exp(beta0)/(1+exp(beta0))
n <- seq(from = 3000, to = 10000, by = 1000)
Power_Proposed1 <- c()
Power_Proposed2 <- c()
for (i in 1:length(n)) {
  Power_Proposed1[i] <- SPCompute:::Compute_Power_Sim(parameters, n = n[i], mode = "dominant", covariate = "continuous", response = "binary")
  Power_Proposed2[i] <- SPCompute:::Compute_Power_Expanded(parameters, n = n[i], mode = "dominant", covariate = "continuous", response = "binary")
}

### First implementation of Demidenko (dichotomize)
set.seed(123)
sim_n <- 300000
G <- sample(c(0,1), size = sim_n, replace = T, prob = c((1-parameters$pG)^2, (1 - ((1-parameters$pG)^2))))
E <- others_para$gamma0 + parameters$gammaG * G + rnorm(sim_n, sd = others_para$sigmaError)
E_new <- ifelse(E > parameters$muE, 1, 0)
y <- beta0 + parameters$betaG * G + parameters$betaE*E + rlogis(sim_n)
y <- ifelse(y>=0, 1, 0)
summary(glm(y~G+E_new, family = binomial()))
summary(glm(E_new~G, family = binomial()))
#### on the misspecified model: betaE = log(3.94), gammaG = log(0.31)

sim1_data_case2_givenSize <- data.frame(n = n)
sim1_data_case2_givenSize$n <- n
sim1_data_case2_givenSize$Demidenko_dicho <- c(0.664, 0.786, 0.868, 0.921, 0.954, 0.973, 0.985, 0.992)
sim1_data_case2_givenSize$Demidenko_noE <- c(0.618, 0.742, 0.831, 0.892, 0.932, 0.958, 0.975, 0.985)
sim1_data_case2_givenSize$Quanto <- c(0.642, 0.765, 0.851, 0.908, 0.944, 0.967, 0.981, 0.989)
sim1_data_case2_givenSize$Proposed1 <- Power_Proposed1
sim1_data_case2_givenSize$Proposed2 <- Power_Proposed2

set.seed(123)
emp_powers <- c()
for (i in 1:length(n)) {
  correct <- c()
  for (j in 1:1000) {
    k <- n[i]
    G <- sample(c(0,1), size = k, replace = T, prob = c((1-parameters$pG)^2, (1 - ((1-parameters$pG)^2))))
    E <- others_para$gamma0 + parameters$gammaG * G + rnorm(n = k, sd = others_para$sigmaError)
    y <- others_para$beta0 + parameters$betaG*G + parameters$betaE*E + rlogis(k)
    y <- ifelse(y>0, 1, 0)
    mod <- glm(y~G+E, family = binomial(link = "logit"))
    correct[j] <- summary(mod)$coefficients[2,4] <= 0.05
  }
  emp_powers[i] <- mean(correct)
}
sim1_data_case2_givenSize$Empirical <- emp_powers
save(file = "sim1_data_case2_givenSize.rda", sim1_data_case2_givenSize)
pdf(file = "sim1_case2_power", height = 4, width = 5)
plot(emp_powers~n, type = 'l', ylab = "Estimated Power", ylim = c(0.5,1), xlab = "Sample Size", lty = "solid")
lines(Power_Proposed1~n, col = "red", lty = "dashed")
lines(Power_Proposed2~n, col = "blue", lty = "dotted")
lines(sim1_data_case2_givenSize$Demidenko_dicho~n, col = "green", lty = "dotdash")
lines(sim1_data_case2_givenSize$Demidenko_noE~n, col = "pink", lty = "dotdash")
lines(sim1_data_case2_givenSize$Quanto~n, col = "purple", lty = "longdash")
legend(x = 6000, y = 0.75,   # Coordinates (x also accepts keywords)
       legend = c("Oracle", "P1.SS", "P2.RD", "Demidenko","Demidenko*", "Quanto"), # Vector with the name of each group
       col = c("black", "red", "blue", "green", "pink", "purple"), bty = "n", cex = 0.8,
       lty = c("solid", "dashed", "dotted", "dotdash", "dotdash", "longdash")
)
dev.off()

methods <- c("Pro1", "Pro2",  "Demidenko", "Demidenko*", "Quanto")
# CR <- c(mean(upper >= Power_Proposed1 & Power_Proposed1 >= lower), mean(upper >= Power_Proposed2 & Power_Proposed2 >= lower), mean(upper >= sim1_data_case1_givenSize$Demidenko & sim1_data_case1_givenSize$Demidenko >= lower), mean(upper >= sim1_data_case1_givenSize$Quanto & sim1_data_case1_givenSize$Quanto >= lower))
MAE <- c(mean(abs(emp_powers - Power_Proposed1)), mean(abs(emp_powers - Power_Proposed2)), mean(abs(emp_powers -sim1_data_case2_givenSize$Demidenko_dicho)), mean(abs(emp_powers -sim1_data_case2_givenSize$Demidenko_noE)), mean(abs(emp_powers -sim1_data_case2_givenSize$Quanto)))
MaxAE <- c(max(abs(emp_powers - Power_Proposed1)), max(abs(emp_powers - Power_Proposed2)), max(abs(emp_powers -sim1_data_case2_givenSize$Demidenko_dicho)), max(abs(emp_powers -sim1_data_case2_givenSize$Demidenko_noE)) ,max(abs(emp_powers -sim1_data_case2_givenSize$Quanto)))
summary_data <- data.frame(methods = methods, MAE = MAE, MaxAE = MaxAE)
summary_data

sim1_data_case2_givenPower <- data.frame(RG = seq(1.1,2.5, by = 0.2))
n_Proposed1 <- c()
n_Proposed2 <- c()
set.seed(123)
Pry <- c()
for (i in 1:length(sim1_data_case2_givenPower$RG)) {
  parameters$betaG <- log(sim1_data_case2_givenPower$RG[i])
  others_para <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "continuous")
  beta0 <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "continuous", seed = 2)$beta0
  Pry[i] <- exp(beta0)/(1+exp(beta0))
  n_Proposed1[i] <- SPCompute:::Compute_Size_Sim(parameters, PowerAim = 0.8, mode = "dominant", covariate = "continuous", response = "binary")
}
for (i in 1:length(sim1_data_case2_givenPower$RG)) {
  parameters$betaG <- log(sim1_data_case2_givenPower$RG[i])
  others_para <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "continuous")
  beta0 <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "continuous", seed = 2)$beta0
  Pry[i] <- exp(beta0)/(1+exp(beta0))
  n_Proposed2[i] <- SPCompute:::Compute_Size_Expanded(parameters, PowerAim = 0.8, mode = "dominant", covariate = "continuous", response = "binary", lower.lim.n = 300)
}
sim1_data_case2_givenPower$Pro1 <- ceiling(n_Proposed1)
sim1_data_case2_givenPower$Pro2 <- ceiling(n_Proposed2)
sim1_data_case2_givenPower$Quanto <- 5 * c(6858, 872, 355, 202, 136, 100, 78, 64)
sim1_data_case2_givenPower$Demidenko_D <- c(33673, 4251, 1726, 1009, 678, 501, 394, 324)
sim1_data_case2_givenPower$Demidenko_NE <- c(39689, 4820, 1892, 1096, 718, 518, 400, 322)
sim1_data_case2_givenPower$oracle <- c()
for (i in 1:length(sim1_data_case2_givenPower$RG)) {
  set.seed(123)
  sim1_data_case2_givenPower$oracle[i] <- compute_oracle_size(OR = sim1_data_case2_givenPower$RG[i], setting = 2)
}
save(file = "sim1_data_case2_givenPower.rda", sim1_data_case2_givenPower)
pdf(file = "sim1_case2_ss.pdf", height = 4, width = 5)
plot(sim1_data_case2_givenPower$Pro1 ~ sim1_data_case2_givenPower$RG, type = 'l', ylab = "Estimated Sample Size", col = "red", xlab = "OR", lty = "dashed")
lines(sim1_data_case2_givenPower$oracle ~ sim1_data_case2_givenPower$RG, col = "black")
lines(sim1_data_case2_givenPower$Pro2 ~ sim1_data_case2_givenPower$RG, col = "blue", lty = "dotted")
lines(sim1_data_case2_givenPower$Demidenko_D ~ sim1_data_case2_givenPower$RG, col = "green", lty = "dotdash")
lines(sim1_data_case2_givenPower$Demidenko_NE ~ sim1_data_case2_givenPower$RG, col = "pink", lty = "dotdash")
lines(sim1_data_case2_givenPower$Quanto ~ sim1_data_case2_givenPower$RG, col = "purple", lty = "longdash")
legend(x = 1.5, y = 40000,   # Coordinates (x also accepts keywords)
       legend = c("Oracle", "P1.SS", "P2.RD","Demidenko", "Demidenko*", "Quanto"), # Vector with the name of each group
       lty = c("solid", "dashed", "dotted", "dotdash", "dotdash", "longdash"),
       col = c("black", "red", "blue" , "green", "pink", "purple"), bty = "n", cex = 0.8
)
dev.off()


############# Comparison number 3: (Case-Control Study, no E) #########################
### We study two cases:
### 1. When sample sizes are given, compare the computed powers with empirical powers
### 2. When powers are given, compute and compare the computed required sample sizes
### Assuming the case to control ratio is 1 to 1 (corresponding to preva = 0.2)
set.seed(123)
parameters <- list(preva = 0.5, pG = 0.1, betaG = log(1.5))
others_para <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "none")
beta0 <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "none", seed = 2)
exp(beta0)/(1+exp(beta0))
N <- seq(from = 300, to = 1400, by = 100)
n <- 2*N
Power_Proposed1 <- c()
Power_Proposed2 <- c()
for (i in 1:length(n)) {
  Power_Proposed1[i] <- SPCompute:::Compute_Power_Sim(parameters, n = n[i], mode = "dominant", covariate = "none", response = "binary")
  Power_Proposed2[i] <- SPCompute:::Compute_Power_Expanded(parameters, n = n[i], mode = "dominant", covariate = "none", response = "binary")
}
sim1_data_case0_givenSize <- data.frame(n = n)
sim1_data_case0_givenSize$n <- n
sim1_data_case0_givenSize$Demidenko <- c(0.487, 0.605, 0.701, 0.778, 0.837, 0.882, 0.916, 0.94, 0.958, 0.971, 0.98, 0.986)
sim1_data_case0_givenSize$Quanto <- c(0.4910, 0.6092, 0.7058, 0.7823, 0.8412, 0.8857, 0.9186, 0.9426, 0.9599, 0.9722, 0.9809, 0.9870)
sim1_data_case0_givenSize$Quanto <- round(sim1_data_case0_givenSize$Quanto, 3)
sim1_data_case0_givenSize$Proposed1 <- round(Power_Proposed1, 3)
sim1_data_case0_givenSize$Proposed2 <- round(Power_Proposed2, 3)
emp_powers <- c()
for (i in 1:length(n)) {
  correct <- c()
  for (j in 1:1000) {
    B <- 10000
    G <- sample(c(0,1), size = B, replace = T, prob = c((1-parameters$pG)^2, (1 - ((1-parameters$pG)^2))))
    y <- beta0 + parameters$betaG*G + rlogis(B)
    y <- ifelse(y>=0, 1, 0)
    data <- data.frame(G, y)
    case_data <- data %>% filter(y == 1) %>% sample_n(N[i], replace = F)
    control_data <- data %>% filter(y == 0) %>% sample_n((n[i] - N[i]), replace = F)
    sampled_data <- rbind(case_data, control_data)
    mod <- glm(y~G, family = binomial(link = "logit"), data = sampled_data)
    correct[j] <- summary(mod)$coefficients[2,4] <= 0.05
  }
  emp_powers[i] <- mean(correct)
}
sim1_data_case0_givenSize$Empirical_CaseControl <- emp_powers
save(file = "sim1_data_case0_givenSize.rda", sim1_data_case0_givenSize)
pdf(file = "sim1_case3_power", height = 4, width = 5)
plot(emp_powers~n, type = 'l', lty = "solid", ylab = "Estimated Power", xlab = "Sample Size", ylim = c(0.5,1))
lines(Power_Proposed1~n, col = "red", lty = "dashed")
lines(Power_Proposed2~n, col = "blue", lty = "dotted")
lines(sim1_data_case0_givenSize$Demidenko~n, col = "green", lty = "dotdash")
lines(sim1_data_case0_givenSize$Quanto~n, col = "purple", lty = "longdash")
legend(x = 1500, y = 0.7,   # Coordinates (x also accepts keywords)
       legend = c("Oracle", "P1.SS", "P2.RD", "Demidenko", "Quanto"), # Vector with the name of each group
       col = c("black", "red", "blue", "green", "purple"),
       lty = c("solid", "dashed", "dotted", "dotdash", "longdash"),
       bty = "n", cex = 0.8
)
dev.off()

### summary of coverage probability, and a summary of error rate;
methods <- c("Pro1", "Pro2", "Demidenko", "Quanto")
# CR <- c(mean(upper >= Power_Proposed1 & Power_Proposed1 >= lower), mean(upper >= Power_Proposed2 & Power_Proposed2 >= lower), mean(upper >= sim1_data_case0_givenSize$Demidenko & sim1_data_case0_givenSize$Demidenko >= lower), mean(upper >= sim1_data_case0_givenSize$Quanto & sim1_data_case0_givenSize$Quanto >= lower))
MAE <- c(mean(abs(emp_powers - Power_Proposed1)), mean(abs(emp_powers - Power_Proposed2)), mean(abs(emp_powers -sim1_data_case0_givenSize$Demidenko)), mean(abs(emp_powers -sim1_data_case0_givenSize$Quanto)))
MaxAE <- c(max(abs(emp_powers - Power_Proposed1)), max(abs(emp_powers - Power_Proposed2)), max(abs(emp_powers -sim1_data_case0_givenSize$Demidenko)), max(abs(emp_powers -sim1_data_case0_givenSize$Quanto)))
summary_data <- data.frame(methods = methods, MAE = MAE, MaxAE = MaxAE)
summary_data

### Compare the three sample size curves (fix the power to be at 0.8, vary the effect size or MAF):
sim1_data_case0_givenPower <- data.frame(RG = seq(1.1,2.5, by = 0.1))
n_Proposed1 <- c()
n_Proposed2 <- c()
set.seed(123)
Pry <- c()
for (i in 1:length(sim1_data_case0_givenPower$RG)) {
  parameters$betaG <- log(sim1_data_case0_givenPower$RG[i])
  others_para <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "none")
  beta0 <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "none", seed = 2)
  Pry[i] <- exp(beta0)/(1+exp(beta0))
  n_Proposed1[i] <- SPCompute:::Compute_Size_Sim(parameters, PowerAim = 0.8, mode = "dominant", covariate = "none", response = "binary")
}
sim1_data_case0_givenPower$Pro1 <- ceiling(n_Proposed1)
for (i in 1:length(sim1_data_case0_givenPower$RG)) {
  parameters$betaG <- log(sim1_data_case0_givenPower$RG[i])
  others_para <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "none")
  beta0 <- convert_preva_to_intercept(parameters, mode = "dominant", covariate = "none", seed = 2)
  Pry[i] <- exp(beta0)/(1+exp(beta0))
  n_Proposed2[i] <- SPCompute:::Compute_Size_Expanded(parameters, PowerAim = 0.8, mode = "dominant", covariate = "none", response = "binary", lower.lim.n = 200)
}
sim1_data_case0_givenPower$Pro1 <- ceiling(n_Proposed1)
sim1_data_case0_givenPower$Pro2 <- ceiling(n_Proposed2)
sim1_data_case0_givenPower$Quanto <- 2 * c(11235, 3075, 1489, 908, 627, 469, 369, 302, 255, 219, 192, 171, 154, 140, 128)
sim1_data_case0_givenPower$Demidenko <- c(22465, 6152, 2985, 1827, 1268, 953, 755, 622, 528, 458, 405, 363, 330, 302, 280)
sim1_data_case0_givenPower$oracle <- c()
for (i in 1:length(sim1_data_case0_givenPower$RG)) {
  set.seed(123)
  sim1_data_case0_givenPower$oracle[i] <- compute_oracle_size(OR = sim1_data_case0_givenPower$RG[i], setting = 3)
}
save(file = "sim1_data_case0_givenPower.rda", sim1_data_case0_givenPower)
pdf(file = "sim1_case3_ss", height = 4, width = 5)
plot(sim1_data_case0_givenPower$Pro1 ~ sim1_data_case0_givenPower$RG, lty = "dashed", type = 'l', ylab = "Estimated Sample Size", col = "red", xlab = "OR")
lines(sim1_data_case0_givenPower$oracle ~ sim1_data_case0_givenPower$RG)
lines(sim1_data_case0_givenPower$Pro2 ~ sim1_data_case0_givenPower$RG, lty = "dotted", col = "blue")
lines(sim1_data_case0_givenPower$Demidenko ~ sim1_data_case0_givenPower$RG, lty = "dotdash", col = "green")
lines(sim1_data_case0_givenPower$Quanto ~ sim1_data_case0_givenPower$RG, lty = "longdash", col = "purple")
legend(x = 1.8, y = 17000,   # Coordinates (x also accepts keywords)
       legend = c("Oracle","P1.SS", "P2.RD", "Demidenko", "Quanto"), # Vector with the name of each group
       lty = c("solid", "dashed", "dotted", "dotdash", "longdash"),
       col = c("black", "red", "blue" , "green", "purple"), bty = "n", cex = 0.8
)
dev.off()




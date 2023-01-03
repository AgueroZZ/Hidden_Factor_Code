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

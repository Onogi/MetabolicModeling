#Analyze simulated data using rrBLUP, quadratic programming, and proposed models

#Simulation serial number
sim <- 1


#Read files########################################################################
#Relationship matrix
G <- as.matrix(read.csv("Gmatrix.csv", header = TRUE))
N <- nrow(G)
colnames(G) <- rownames(G) <- 1:N

#Metabolic network
K <- 350
J <- 336
Bound <- read.csv("Bound.nonzero.csv", header = TRUE, row.names = 1)
Positiveflux <- Bound[1, ] >= 0

Pos_Jc <- which(Positiveflux)
Pos_Jn <- c(1:J)[-c(Pos_Jc)]
Jc <- length(Pos_Jc)
Jn <- length(Pos_Jn)

#Folds of CV
FoldID <- as.matrix(read.csv("netGS/foldid.csv", header = FALSE))

#Simulated data
Y <- as.matrix(read.csv(paste("Simulation", sim, ".Y.csv", sep = ""), header = TRUE, row.names = 1))
M <- as.matrix(read.csv(paste("Simulation", sim, ".M.csv", sep = ""), header = TRUE, row.names = 1))
V <- as.matrix(read.csv(paste("Simulation", sim, ".V.csv", sep = ""), header = TRUE, row.names = 1))


#Cross validation###############################################################
library(rstan)
library(rrBLUP)
library(snow)
library(CVXR)
library(MegaLMM)

#Training and test data
Nfold <- 3
Fold <- matrix(0, nrow = max(table(FoldID[, sim])), Nfold)
for(fold in 1:Nfold){
  Fold[1:table(FoldID[, sim])[fold], fold] <- which(FoldID[, sim] == fold)
}

#rrBLUP
Predict.rrBLUP <- numeric(N)
for(fold in 1:Nfold){
  
  Test <- Fold[, fold]
  Train <- c(1:N)[-Test]
  
  Result.temp <- mixed.solve(Y[Train], diag(N)[Train, ], G)
  Predict.rrBLUP[Test] <- Result.temp$u[Test] + as.numeric(Result.temp$beta)
}


#Quadratic programming
Flux.Bio <- matrix(NA, N, J)
for(i in 1:N){
  
  b <- Variable(J)
  W <- ((i - 1) * K + 1):(i * K)
  Amat <- M[W,]
  Result.temp <- psolve(Problem(Minimize(norm2(Y[i] - b[J])),
                              list(Amat%*%b == 0, 
                                   b[Pos_Jc] >= 0))
  )
  
  Flux.Bio[i, ] <- as.numeric(Result.temp$getValue(b))
}
write.csv(Flux.Bio,paste("Simulation", sim, ".Flux.csv", sep = ""))
#=>These estimates are used for the proposed model

#Predict for new genotypes (non steady state)
Nonss.Bio <- matrix(NA, N, J)
for(fold in 1:Nfold){
  
  Test <- Fold[,fold]
  Train <- c(1:N)[-Test]
  
  for(j in 1:J){
    Result.temp <- mixed.solve(Flux.Bio[Train, j], diag(N)[Train, ], G)
    Nonss.Bio[Test, j] <- Result.temp$u[Test] + as.numeric(Result.temp$beta)
  }
}

#Make fluxes steady state
SS.Bio.Bp1 <- matrix(NA, N, J)
for(i in 1:N){
  p <- Nonss.Bio[i,]
  b <- Variable(J)
  
  W <- ((i - 1) * K + 1):(i * K)
  Amat <- M[W, ]
  
  try(Result.temp <- psolve(Problem(Minimize(norm2(b/p - 1)),
                                    list(Amat%*%b == 0, 
                                         b[Pos_Jc] >= 0))
  ))
  
  if(Result.temp$status == "solver_error"){
    SS.Bio.Bp1[i, ] <- NA
  }else{
    SS.Bio.Bp1[i, ] <- as.numeric(Result.temp$getValue(b))
  }
}
Predict.QP <- SS.Bio.Bp1[, J]


#MegaLMM
n_iter <- 100
Nblock <- 100

#Objects storing results
Predict.Mega <- matrix(0, 67, 3)
colnames(Predict.Mega) <- c(50, 80, 160)
MegaLMM_state <- as.list(numeric(3))
names(MegaLMM_state) <- c(50, 80, 160)
for(k in c(50, 80, 160)){
  MegaLMM_state[[as.character(k)]] <- as.list(numeric(3))#Number of folds
}

#Read files
Flux <- as.matrix(read.csv(paste("Simulation", sim, ".Flux.csv", sep = ""), header = TRUE, row.names = 1))

for(k in c(50, 80, 160)){
  run_parameters <- MegaLMM_control(
    max_NA_groups = 3,
    scale_Y = TRUE,
    h2_divisions = 20,
    h2_step_size = NULL,
    burn = 0,
    K = k
  )
  
  priors <- MegaLMM_priors(
    tot_Y_var = list(V = 0.5,   nu = 5),
    tot_F_var = list(V = 18/20, nu = 20),
    Lambda_prior = list(
      sampler = sample_Lambda_prec_horseshoe,
      prop_0 = 0.1,
      delta = list(shape = 3, scale = 1),
      delta_iterations_factor = 100
    ),
    h2_priors_resids_fun = function(h2s,n) 1,
    h2_priors_factors_fun = function(h2s,n) 1
  )
  
  #Run
  for(fold in 1:Nfold){
    
    Test <- sort(Fold[, fold])
    Train <- c(1:N)[-Test]
    m <- mean(Y[Train])
    s <- sd(Y[Train])
    
    Resvariable <- cbind(Flux, Y)
    colnames(Resvariable) <- 1:ncol(Resvariable)
    rownames(Resvariable) <- 1:nrow(Resvariable)
    Resvariable[Test, ] <- NA
    data <- data.frame(Intercept = 1, ID = as.character(1:nrow(Resvariable)))
    
    MegaLMM_state[[as.character(k)]][[fold]] <- setup_model_MegaLMM(Resvariable,
                                                                    ~(1|ID),
                                                                    data = data,
                                                                    relmat = list(ID = G),
                                                                    run_parameters=run_parameters,
                                                                    run_ID = paste("K", k, "_Simulation", sim, "_fold", fold, sep = "")
    )
    maps <- make_Missing_data_map(MegaLMM_state[[as.character(k)]][[fold]])
    MegaLMM_state[[as.character(k)]][[fold]] <- set_Missing_data_map(MegaLMM_state[[as.character(k)]][[fold]], maps$Missing_data_map)
    MegaLMM_state[[as.character(k)]][[fold]] <- set_priors_MegaLMM(MegaLMM_state[[as.character(k)]][[fold]], priors)
    MegaLMM_state[[as.character(k)]][[fold]] <- initialize_variables_MegaLMM(MegaLMM_state[[as.character(k)]][[fold]])
    MegaLMM_state[[as.character(k)]][[fold]] <- initialize_MegaLMM(MegaLMM_state[[as.character(k)]][[fold]])
    MegaLMM_state[[as.character(k)]][[fold]] <- clear_Posterior(MegaLMM_state[[as.character(k)]][[fold]])
    
    for(i  in 1:Nblock) {
      print(sprintf('Run %d',i))
      MegaLMM_state[[as.character(k)]][[fold]] <- sample_MegaLMM(MegaLMM_state[[as.character(k)]][[fold]], n_iter)
      MegaLMM_state[[as.character(k)]][[fold]] <- save_posterior_chunk(MegaLMM_state[[as.character(k)]][[fold]])
      print(MegaLMM_state[[as.character(k)]][[fold]])
      if(i < 91) {
        MegaLMM_state[[as.character(k)]][[fold]] <- reorder_factors(MegaLMM_state[[as.character(k)]][[fold]], drop_cor_threshold = 0.6)
        MegaLMM_state[[as.character(k)]][[fold]] <- clear_Posterior(MegaLMM_state[[as.character(k)]][[fold]])
        print(MegaLMM_state[[as.character(k)]][[fold]]$run_parameters$burn)
      }
    }
    plot(MegaLMM_state[[as.character(k)]][[fold]])
    
    MegaLMM_state[[as.character(k)]][[fold]]$Posterior <- reload_Posterior(MegaLMM_state[[as.character(k)]][[fold]])
    U_hat <- numeric(N)
    for(ite in 1:500){
      Intercept <- MegaLMM_state[[as.character(k)]][[fold]]$Posterior$B1[ite, , ]
      Intercept <- matrix(Intercept, N, J + 1, byrow=TRUE)
      Predict.temp <- Intercept + MegaLMM_state[[as.character(k)]][[fold]]$Posterior$U_F[ite, , ] %*% 
        MegaLMM_state[[as.character(k)]][[fold]]$Posterior$Lambda[ite, , ] + 
        MegaLMM_state[[as.character(k)]][[fold]]$Posterior$U_R[ite, , ]
      
      Predict.temp <- Predict.temp * s + m
      
      U_hat <- U_hat + Predict.temp[, J + 1]
    }
    
    Predict.Mega[Test, as.character(k)] <- as.numeric(U_hat[Test])/500
  }#fold
}#k


#Proposed model
Weight <- c(0.16, 0.32, 0.48, 0.64, 0.80, 0.96)
Model10rho <- stan_model(file = "MetabolicModeling10rho.stan")

#Create Data object
Predict.MM <- matrix(NA, N, length(Weight))
colnames(Predict.MM) <- Weight
Seed <- sample.int(.Machine$integer.max, 1)

Data <- as.list(numeric(Nfold))
for(fold in 1:Nfold){
  
  Test <- sort(Fold[,fold])
  Train <-c(1:N)[-Test]
  
  Nn <- length(Train)
  Pos_NonMiss <- Train
  
  m <- mean(Y[Train])
  s <- sd(Y[Train])
  Y.scale <- as.vector((Y - m)/s)
  
  #Calculate variance components of fluxes
  Variance.Train <- matrix(0, nrow = J, ncol = 3)
  colnames(Variance.Train) <- c("Ave", "Vu", "Ve")
  for(j in 1:J){
    Y.flux <- as.numeric(Flux[Train, j])
    temp <- mixed.solve(y = Y.flux, K = G[Train, Train])
    Variance.Train[j, ] <- c(temp$beta,
                             temp$Vu,
                             temp$Ve)
  }
  Variance.Train <- data.frame(Variance.Train)
  
  #Scale the variance of V
  Delta <- sqrt(Variance.Train$Vu + Variance.Train$Ve)
  M.delta <-t(t(M) * Delta)
  
  #Create Alpha and Ma
  Alpha <- 10 - Variance.Train$Ave/Delta
  Ma <- matrix(0, N, K)
  for(i in 1:N){
    for(k in 1:K){
      W <-(i - 1) * K + k
      Ma[i, k] <- sum(M.delta[W, ] * Alpha)
    }
  }
  
  #Creat Muv
  Muv <- Variance.Train$Ave/Delta
  
  Data[[fold]] <- list(J = J, Jc = Jc, Jn = Jn,
                       K = K, N = N, Nn = Nn,
                       M = M.delta, G = G, Y = Y.scale,
                       W = NULL,
                       Muv = Muv,
                       Alpha = Alpha,
                       Ma = Ma,
                       Pos_Jc = Pos_Jc, Pos_Jn = Pos_Jn,
                       Pos_NonMiss = Train)
}

#Add noise
Noise <- as.list(numeric(Nfold))
for(fold in 1:Nfold){
  Noise[[fold]] <- matrix(rnorm(N * K, 0, sqrt(0.0001)), N, K)
}
for(fold in 1:Nfold){
  Data[[fold]]$Ma <- Data[[fold]]$Ma + Noise[[fold]]
}

#Run
Result.MM <- as.list(numeric(length(Weight)))
names(Result.MM) <- Weight
for(W in Weight){
  
  cat(W, "\n")
  
  cl <- makeCluster(Nfold, type = "SOCK")
  clusterExport(cl, c("Model10rho",
                      "sampling",
                      "summary",
                      "Seed"))
  
  for(fold in 1:Nfold) Data[[fold]]$W <- W
  
  Result.MM[[as.character(W)]] <- as.list(numeric(Nfold))
  Result.MM[[as.character(W)]] <- parLapply(cl, 
                                         Data, 
                                         function(Data){
                                           summary(
                                             sampling(
                                               object = Model10rho,
                                               data = Data,
                                               iter = 5000,
                                               warmup = 4000,
                                               thin = 1,
                                               chains = 1,
                                               algorithm = "NUTS",
                                               seed = Seed
                                             )
                                           )
                                         }
  )
  stopCluster(cl)
  save.image(paste("Simulation", sim, ".RData", sep = ""))
}

#Extract results and predict
Predict.MM <- matrix(NA, N, length(Weight))
colnames(Predict.MM) <- Weight
for(W in Weight){
  
  cat(W,"\n")
  for(fold in 1:Nfold){
    
    Test <- Fold[, fold]
    Train <- c(1:N)[-Test]
    
    m <- mean(Y[Train])
    s <- sd(Y[Train])
    
    cat(fold, "\n")
    
    #Estimates of V
    Result.MM.Vhat <- matrix(0, nrow = N, ncol = J)
    start <- 2*J+2*J*N+6+J*N
    for(j in 1:J){
      Result.MM.Vhat[, j] <- Result.MM[[as.character(W)]][[fold]]$summary[((j - 1) * N + start):(j * N + start - 1), "mean"]
    }
    
    #Intercept and slope
    Ahat <- Result.MM[[as.character(W)]][[fold]]$summary["A", "mean"]
    Bhat <- Result.MM[[as.character(W)]][[fold]]$summary["B", "mean"]
    
    #Predict
    Predict.MM[Test, as.character(W)] <- (Ahat + Bhat * (Result.MM.Vhat[Test, J])) * s + m
  }
}


#Output###############################################################################
#Predicted values
Predict <- data.frame(Y = Y,
                      Fold = FoldID[, sim],
                      rrBLUP = Predict.rrBLUP,
                      QP = Predict.QP,
                      MegaLMM050 = Predict.Mega[, 1],
                      MegaLMM080 = Predict.Mega[, 2],
                      MegaLMM160 = Predict.Mega[, 3],
                      MM016 = Predict.MM[,1],
                      MM032 = Predict.MM[,2],
                      MM048 = Predict.MM[,3],
                      MM064 = Predict.MM[,4],
                      MM080 = Predict.MM[,5],
                      MM096 = Predict.MM[,6])
write.csv(Predict, paste("Simulation", sim, "_Predict.csv", sep = ""), row.names = FALSE)


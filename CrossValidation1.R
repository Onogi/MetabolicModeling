#Conduct cross validation for the real data
#rrBLUP, quadratic programming, and proposed models are compared.

#CV number
CV <- 1


#Read files######################################################################
#Relationship matrix
G <- as.matrix(read.csv("Gmatrix.csv", header = TRUE))
N <- nrow(G)
colnames(G) <- rownames(G) <- 1:N

#Stoichiometry matrix
Smatrix <- read.csv("Smatrix.nonzero.csv", header = TRUE, row.names = 1)
dim(Smatrix)
#350 336
K <- nrow(Smatrix)
J <- ncol(Smatrix)

Bound <- read.csv("Bound.nonzero.csv", header = TRUE, row.names = 1)
Positiveflux <- Bound[1, ] >= 0

#Constraints
Pos_Jc <- which(Positiveflux)
Pos_Jn <- c(1:J)[-c(Pos_Jc)]
Jc <- length(Pos_Jc)
Jn <- length(Pos_Jn)

#Variety-unique Jth flux (biomass reaction)
MJ <- read.csv("Biomass.nonzero.csv", header = TRUE, row.names = 1)

#Estimated flux for each variety
Flux <- read.csv("TongH2020SupplementaryData2.csv", header = TRUE, stringsAsFactors = FALSE)

#Observed biomass
Y <- as.numeric(unlist(read.csv("netGS/biomass_optN.csv", header = FALSE, row.names = 1)))

#Folds of CV
FoldID <- as.matrix(read.csv("netGS/foldid.csv", header = FALSE))


#Variance component estimates##################################################
library(rrBLUP)
Variance <- matrix(0, nrow = J, ncol = 6)
colnames(Variance) <- c("Ave", "CV", "Vp", "Vu", "Ve", "h2")
for(j in 1:J){
  y <- as.numeric(Flux[j, -c(1:4)])
  temp <- mixed.solve(y = y, K = G)
  Variance[j, ] <- c(temp$beta,
                     abs(sd(y)/temp$beta),
                     var(y),
                     temp$Vu,
                     temp$Ve,
                     temp$Vu/(temp$Vu + temp$Ve))
}
Variance <- data.frame(Variance)
rm(y, temp, j)


#Do CV###########################################################################
library(MegaLMM)
library(CVXR)
#=>Matlab version of CVX was used in Tong et al. (2020)

#Training and test data
Nfold <- 3
Fold <- matrix(0, nrow = max(table(FoldID[, CV])), Nfold)
for(fold in 1:Nfold){
  Fold[1:table(FoldID[, CV])[fold], fold] <- which(FoldID[, CV] == fold)
}

#rrBLUP
Predict.rrBLUP <- numeric(N)
for(fold in 1:Nfold){
  Test <- sort(Fold[,fold])
  Train <-c(1:N)[-Test]
  
  m <- mean(Y[Train])
  s <- sd(Y[Train])
  Y.scale <- (Y - m)/s
  
  Result.rrBLUP <- mixed.solve(Y.scale[Train], diag(N)[Train, ], G)
  Predict.rrBLUP[Test] <- (Result.rrBLUP$u[Test] + as.numeric(Result.rrBLUP$beta)) * s + m
}

#Quadratic programming
Predict.QP.nonss <- matrix(NA, N, J)#non steady state
for(fold in 1:Nfold){
  
  Test <- Fold[,fold]
  Train <-c(1:N)[-Test]
  
  for(j in 1:J){
    Result.temp <- mixed.solve(as.numeric(Flux[j, -c(1:4)][Train]), diag(N)[Train, ], G)
    Predict.QP.nonss[Test, j] <- Result.temp$u[Test] + as.numeric(Result.temp$beta)
  }
}

Predict.QP <- matrix(NA, N, J)#steady state
for(i in 1:N){
  p <- Predict.QP.nonss[i, ]
  b <- Variable(J)
  
  Amat <- as.matrix(Smatrix)
  Amat[, J] <- MJ[, i]
  
  Result.temp <- psolve(Problem(Minimize(norm2(b/p - 1)),
                                list(Amat%*%b == 0, 
                                     b[Pos_Jc] >= 0))
  )
  
  Predict.QP[i, ] <- as.numeric(Result.temp$getValue(b))
}


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
    
    Resvariable <- cbind(t(Flux[, -c(1:4)]), Y)
    colnames(Resvariable) <- 1:ncol(Resvariable)
    rownames(Resvariable) <- 1:nrow(Resvariable)
    Resvariable[Test, ] <- NA
    data <- data.frame(Intercept = 1, ID = as.character(1:nrow(Resvariable)))
    #=>Intercept is unnecessary
    
    MegaLMM_state[[as.character(k)]][[fold]] <- setup_model_MegaLMM(Resvariable,
                                                       ~(1|ID),
                                                       data = data,
                                                       relmat = list(ID = G),
                                                       run_parameters=run_parameters,
                                                       run_ID = paste("K", k, "_CV", CV, "_fold", fold, sep = "")
    )
    maps <- make_Missing_data_map(MegaLMM_state[[as.character(k)]][[fold]])
    MegaLMM_state[[as.character(k)]][[fold]] <- set_Missing_data_map(MegaLMM_state[[as.character(k)]][[fold]], maps$Missing_data_map)
    MegaLMM_state[[as.character(k)]][[fold]] <- set_priors_MegaLMM(MegaLMM_state[[as.character(k)]][[fold]], priors)
    MegaLMM_state[[as.character(k)]][[fold]] <- initialize_variables_MegaLMM(MegaLMM_state[[as.character(k)]][[fold]])
    MegaLMM_state[[as.character(k)]][[fold]] <- initialize_MegaLMM(MegaLMM_state[[as.character(k)]][[fold]])
    MegaLMM_state[[as.character(k)]][[fold]] <- clear_Posterior(MegaLMM_state[[as.character(k)]][[fold]])
    
    for(i  in 1:Nblock) {
      print(sprintf('Run %d', i))
      MegaLMM_state[[as.character(k)]][[fold]] <- sample_MegaLMM(MegaLMM_state[[as.character(k)]][[fold]], n_iter)
      MegaLMM_state[[as.character(k)]][[fold]] <- save_posterior_chunk(MegaLMM_state[[as.character(k)]][[fold]])
      print(MegaLMM_state[[as.character(k)]][[fold]])
      if(i < 91) {
        MegaLMM_state[[as.character(k)]][[fold]] <- reorder_factors(MegaLMM_state[[as.character(k)]][[fold]], drop_cor_threshold = 0.6)
        MegaLMM_state[[as.character(k)]][[fold]] <- clear_Posterior(MegaLMM_state[[as.character(k)]][[fold]])
        print(MegaLMM_state[[as.character(k)]][[fold]]$run_parameters$burn)
      }
    }
    #plot(MegaLMM_state[[as.character(k)]][[fold]])
    
    MegaLMM_state[[as.character(k)]][[fold]]$Posterior <- reload_Posterior(MegaLMM_state[[as.character(k)]][[fold]])
    U_hat <- numeric(N)
    for(ite in 1:500){
      Intercept <- MegaLMM_state[[as.character(k)]][[fold]]$Posterior$B1[ite, , ]
      Intercept <- matrix(Intercept, N, J + 1, byrow = TRUE)
      Predict.temp <- Intercept + MegaLMM_state[[as.character(k)]][[fold]]$Posterior$U_F[ite, , ] %*% 
        MegaLMM_state[[as.character(k)]][[fold]]$Posterior$Lambda[ite, , ] + 
        MegaLMM_state[[as.character(k)]][[fold]]$Posterior$U_R[ite, , ]
      
      Predict.temp <- Predict.temp * s + m
      
      U_hat <- U_hat + Predict.temp[, J + 1]
    }
    
    Predict.Mega[Test, as.character(k)] <- as.numeric(U_hat[Test])/500
  }#fold
}#k


#Proposed model based on metabolic networks
Weight<-c(0.16, 0.32, 0.48, 0.64, 0.80, 0.96)
library(rstan)
library(rrBLUP)
library(snow)

#Use Model10rho
#change prior
Model10rho <- stan_model(file = "MetabolicModeling10rho.stan")

#Metabolite modeling
Predict.MM <- matrix(NA, N, length(Weight))
colnames(Predict.MM) <- Weight
Seed<-sample.int(.Machine$integer.max, 1)

Data <- as.list(numeric(Nfold))
for(fold in 1:Nfold){
  
  Test <- sort(Fold[,fold])
  Train <-c(1:N)[-Test]
  
  Nn <- length(Train)
  Pos_NonMiss <- Train
  
  m <- mean(Y[Train])
  s <- sd(Y[Train])
  Y.scale <- (Y - m)/s
  
  #Calculate variance components of fluxes
  Variance.Train <- matrix(0, nrow = J, ncol = 3)
  colnames(Variance.Train) <- c("Ave", "Vu", "Ve")
  for(j in 1:J){
    Y.flux <- as.numeric(Flux[j, -c(1:4)])
    Y.flux <- Y.flux[Train]
    temp <- mixed.solve(y = Y.flux, K = G[Train, Train])
    Variance.Train[j, ] <- c(temp$beta,
                             temp$Vu,
                             temp$Ve)
  }
  Variance.Train <- data.frame(Variance.Train)
  
  #Scale the variance of V
  M<-NULL
  for(i in 1:N){
    temp <- Smatrix
    temp[, J] <- MJ[, i]
    M <- rbind(M, temp)
  }
  dim(M)
  #23450   336
  Delta <- sqrt(Variance.Train$Vu + Variance.Train$Ve)
  M <- t(t(M) * Delta)
  dim(M)
  #23450   336
  
  #Create Alpha and Ma
  Alpha <- 10 - Variance.Train$Ave/Delta
  Ma <- matrix(0, N, K)
  for(i in 1:N){
    for(k in 1:K){
      W <- (i - 1) * K + k
      Ma[i, k] <- sum(M[W, ] * Alpha)
    }
  }
  
  #Creat Muv
  Muv <- Variance.Train$Ave/Delta
  
  Data[[fold]] <- list(J = J, Jc = Jc, Jn = Jn,
                     K = K, N = N, Nn = Nn,
                     M = M, G = G, Y = Y.scale,
                     W = NULL,
                     Muv = Muv,
                     Alpha = Alpha,
                     Ma = Ma,
                     Pos_Jc = Pos_Jc, Pos_Jn = Pos_Jn,
                     Pos_NonMiss = Train)
}
#remove objects unique to fold
rm(Nn, M, Y.scale, Muv, Alpha, Ma)
save.image(paste("CrossValidation", CV, ".RData", sep = ""))

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
  save.image(paste("CrossValidation", CV, ".RData", sep = ""))
}

#Extract results and predict
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
                      Fold = FoldID[, CV],
                      rrBLUP = Predict.rrBLUP,
                      QP = Predict.QP[, J],
                      MegaLMM050 = Predict.Mega[, 1],
                      MegaLMM080 = Predict.Mega[, 2],
                      MegaLMM160 = Predict.Mega[, 3],
                      MM016 = Predict.MM[, 1],
                      MM032 = Predict.MM[, 2],
                      MM048 = Predict.MM[, 3],
                      MM064 = Predict.MM[, 4],
                      MM080 = Predict.MM[, 5],
                      MM096 = Predict.MM[, 6])
write.csv(Predict, paste("CV", CV, "_Predict.csv", sep = ""), row.names = FALSE)

#Status of filling constraints
FittingToMa <- NULL
for(W in c(0.16, 0.32, 0.48, 0.64, 0.80, 0.96)){
  
  Ma.x <- NULL
  Ma.y <- NULL
  Y.x <- NULL
  Y.y <- NULL
  
  Result.temp <- Result.MM[[as.character(W)]]
  
  for(fold in 1:Nfold){
    
    Test <- Fold[,fold]
    Train <- c(1:N)[-Test]
    
    x <- as.numeric(t(Data[[fold]]$Ma))
    y <- Result.temp[[fold]]$summary[grep("Mv", rownames(Result.temp[[fold]]$summary)), "mean"]
    Ma.x <- c(Ma.x, x)
    Ma.y <- c(Ma.y, y)
  }
  FittingToMa <- cbind(FittingToMa, cbind(Ma.x, Ma.y))
}

colnames(FittingToMa) <- paste(rep(c("Ma", "Fit"), 6), 
                               rep(c(0.16, 0.32, 0.48, 0.64, 0.80, 0.96), each = 2), sep = "_")
FittingToMa <- data.frame(FittingToMa)
write.csv(FittingToMa, paste("CV", CV, "_FittingToMa.csv", sep = ""), row.names = FALSE)


#Simulated data following the real data

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


#Simulate data based on the real data###########################################
library(MASS)
Smatrix <- as.matrix(Smatrix)
for(sim in 1:2){
  
  #For metabolites of which M[,J] are 0, the last non-zero column is used for adjustment
  U <- matrix(0, N, J)
  V <- matrix(0, N, J)
  
  #Use the estimated variance components
  for(j in 1:J){
    if(Positiveflux[j]){
      repeat{
        U[, j] <- mvrnorm(1, rep(Variance$Ave[j], N), G * Variance$Vu[j])
        V[, j] <- U[, j] + rnorm(N, 0, sqrt(Variance$Ve[j]))
        if(sum(V[, j] > 0) == N) break
      }
    }else{
      U[, j] <- mvrnorm(1, rep(Variance$Ave[j], N), G * Variance$Vu[j])
      V[, j] <- U[, j] + rnorm(N, 0, sqrt(Variance$Ve[j]))
    }
  }
  
  #For each k ,the last column is overwritten
  M <- matrix(0, N * K, J)
  for(i in 1:N){
    for(k in 1:K){
      Last <- max(which(Smatrix[k, ] != 0))
      W <- (i - 1) * K + k
      
      M[W, -Last] <- Smatrix[k, -Last]
      
      v <- M[W, -Last] %*% matrix(V[i, -Last], ncol = 1)
      M[W, Last] <- (-v) / V[i, Last]
    }
  }
  
  A <- 0
  B <- 1
  Y <- A + B * V[, J]
  E <- rnorm(N, 0, sqrt(var(Y)*0.25))
  Y <- Y + E
  
  write.csv(U, paste("Simulation", sim, ".U.csv", sep = ""))
  write.csv(V, paste("Simulation", sim, ".V.csv", sep = ""))
  write.csv(M, paste("Simulation", sim, ".M.csv", sep = ""))
  write.csv(Y, paste("Simulation", sim, ".Y.csv", sep = ""))
}#sim




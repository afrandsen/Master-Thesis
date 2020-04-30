coef

Phi0_temp <- coef[ , ncol(coef)]
Phi1_temp <- coef[ , 1:(ncol(coef)-1)]

z <- data_sub[ , 1:(ncol(data_sub))]
z <- as.matrix(z)
z2 <- t(z)

port_calc <- function(gamma = 5, K = 100, n= 4){
m <- nrow(Phi1_temp)

H <- diag(x = 1, m, m)

H1 <- H[1,]
H1 <- t(as.matrix(H1))

Hx <- H[2:n, ]
Hx <- as.matrix(Hx)

Phi0 <- Phi0_temp
Phi0 <- as.matrix(Phi0)

Phi1 <- Phi1_temp
Phi1 <- as.matrix(Phi1)

Phi0X <- Hx %*% Phi0
Phi1X <- Hx %*% Phi1
Phi01 <- H1 %*% Phi0
Phi11 <- H1 %*% Phi1

Sigma <- cov
SigmaXX <- Hx %*% Sigma %*% t(Hx)
SigmaX2 <- diag(SigmaXX)
SigmaX <- Hx %*% Sigma
SigmaX1 <- Hx %*% Sigma %*% t(H1)
Sigma11 <- H1 %*% Sigma %*% t(H1)

A0 <- matrix(0, nrow = n-1, ncol = 1)
A1 <- matrix(0, nrow = n-1, ncol = m)
B1 <- matrix(0, nrow = 1, ncol = m)
B2 <- matrix(0, nrow = m, ncol = m)
Lambda <- matrix(0, nrow = m, ncol = m)
Gamma <- matrix(0, nrow = m, ncol = m)
Xi <- matrix(0, nrow = m, ncol = m)

for (i in 1:(K-1)) {
  if (i==1) {
    # A0(1).
    A0 <- (1/gamma) * solve(SigmaXX) %*% (Phi0X + 1/2 * SigmaX2 + (1 - gamma) * SigmaX1)
    
    # A1(1).
    A1 <- (1/gamma) * solve(SigmaXX) %*% Phi1X
    
    # B1(1).
    B1 <- Phi11 + t(A0) %*% (Phi1X - gamma * SigmaXX %*% A1) + t(Phi0X + 1/2 * SigmaX2 + (1 - gamma) * SigmaX1) %*% A1
    
    # B2(1).
    B2 <- t(A1) %*% (Phi1X - gamma/2 * SigmaXX %*% A1)
  } else {
    
    # Gem sidste periode A0 og A1.
    A0a <- A0
    A1a <- A1
    
    # Omega.
    Omega <- solve(solve(Sigma)-2 * (1 - gamma) * B2)
    OmegaXX <- Hx %*% Omega %*% t(Hx)
    OmegaX1 <- Hx %*% Omega %*% t(H1)
    OmegaX <- Hx %*% Omega
    
    # A0.
    A0 <- solve(SigmaXX - (1 - gamma) * OmegaXX) %*% (Phi0X + 1/2 * SigmaX2 + (1 - gamma) * (OmegaX1 + OmegaX %*% (t(B1) + 2 * B2 %*% Phi0)))
    
    # A1.
    A1 <- solve(SigmaXX - (1 - gamma) * OmegaXX) %*% (Phi1X + 2 * (1 - gamma) * OmegaX %*% B2 %*% Phi1)
    
    # Lambda.
    Lambda <- B2 %*% Omega %*% t(B2)
    
    # Gamma.
    Gamma <- 2 * B2 %*% Omega
    
    # Xi.
    Xi <- t(Gamma)
    Xix <- Hx %*% Xi
    Xi1 <- H1 %*% Xi
    
    # Beregn B1 og B2 for næste periode.
    
    # B1.
    B1 <- Phi11 + t(A0a) %*% (Phi1X - (SigmaXX - (1 - gamma) * OmegaXX) %*% A1a) + t(Phi0X + 1/2 * SigmaX2 + (1 - gamma) * OmegaX1) %*% A1a + (B1 + 2 * t(Phi0) %*% B2) %*% Phi1 + (1 - gamma) * (B1 %*% t(OmegaX) + t(Phi0) %*% t(Xix)) %*% A1a + (1 - gamma) * (4 * t(Phi0) %*% Lambda + Xi1 + t(A0a) %*% Xix + B1 %*% Xi) %*% Phi1
    
    # B2.
    B2 <- t(A1a) %*% (Phi1X - gamma/2 * SigmaXX %*% A1a) + t(Phi1) %*% (B2 + 2 * (1 - gamma) * Lambda) %*% Phi1 + (1 - gamma) * t(Phi1) %*% t(Xix) %*% A1a
  }
}

if (K > 1) {
  # A0.
  A0 <- solve(SigmaXX - (1 - gamma) * OmegaXX) %*% (Phi0X + 1/2 * SigmaX2 + (1 - gamma) * (OmegaX1 + OmegaX %*% (t(B1) + 2 * B2 %*% Phi0)))
  
  # A1.
  A1 <- solve(SigmaXX - (1 - gamma) * OmegaXX) %*% (Phi1X + 2 * (1 - gamma) * OmegaX %*% B2 %*% Phi1)
  
} else {
  # A0.
  A0 <- 1/gamma * solve(SigmaXX) %*% (Phi0X + 1/2 * SigmaX2 + (1 - gamma) * SigmaX1)
  
  # A1.
  A1 <- 1/gamma * solve(SigmaXX) %*% Phi1X
}

alphafixedT <- matrix(0, nrow = n-1, ncol = ncol(z2))

for (t in (2:ncol(z2))) {
  alphafixedT[, t] <- A0 + A1 %*% z2[, t-1]
}

MyopicD <- matrix(0, nrow = n-1, ncol = ncol(z2)) 

for (K in (2:ncol(z2))) {
  MyopicD[, K] <- (1/gamma) * solve(SigmaXX) %*% (Hx %*% (Phi0 + Phi1 %*% z2[, K-1]) + 1/2 * SigmaX2 + (1-gamma) * SigmaX1)
}

HedgingD <- alphafixedT - MyopicD

TangencyP <- matrix(0, nrow = n-1, ncol = ncol(z2))

for (t in (2:ncol(z2))){
  TangencyP[, t] <- solve(SigmaXX) %*% Hx %*% (Phi0 + Phi1 %*% z2[, t-1]) + 1/2 * SigmaX2
}

GMV <- -solve(SigmaXX)%*%SigmaX1

alphafixedT <- rbind(1-colSums(alphafixedT), alphafixedT)

res_alphafixedT <- matrix(0, nrow = n, ncol = ncol(z2))

for (t in (2:ncol(z2))) {
  res_alphafixedT[, t] <- c(max(0, alphafixedT[1, t])/sum(max(0, alphafixedT[1,t]),
                                        max(0, alphafixedT[2,t]),
                                        max(0, alphafixedT[3,t]),
                                        max(0, alphafixedT[4,t])),
                            
                            max(0, alphafixedT[2, t])/sum(max(0, alphafixedT[1,t]),
                                                 max(0, alphafixedT[2,t]),
                                                 max(0, alphafixedT[3,t]),
                                                 max(0, alphafixedT[4,t])),
                            
                            max(0, alphafixedT[3, t])/sum(max(0, alphafixedT[1,t]),
                                                 max(0, alphafixedT[2,t]),
                                                 max(0, alphafixedT[3,t]),
                                                 max(0, alphafixedT[4,t])),
                            
                            max(0, alphafixedT[4, t])/sum(max(0, alphafixedT[1,t]),
                                                          max(0, alphafixedT[2,t]),
                                                          max(0, alphafixedT[3,t]),
                                                          max(0, alphafixedT[4,t]))
                          )
}

return(list(alphafixedT=alphafixedT,
            MyopicD=MyopicD,
            HedgingD=HedgingD,
            TangencyP=TangencyP,
            GMV=GMV,
            res_alphafixedT=res_alphafixedT)
       )
}

dynport <- port_calc(gamma = 5, K = 100, n = 4)

KOMP <- data.frame(Tan=c(1-sum(rowMeans(dynport$TangencyP)),
                         rowMeans(dynport$TangencyP)),
                   GMV=c(1-sum(dynport$GMV),
                         dynport$GMV),
                   Myopic=c(1-sum(rowMeans(dynport$MyopicD)),
                            rowMeans(dynport$TangencyP)),
                   IHD=c(1-sum(rowMeans(dynport$HedgingD)), rowMeans(dynport$HedgingD))
                   )

BIG_TABLE <- matrix(0, nrow = 16, ncol = 8)
BIG_R_TABLE <- matrix(0, nrow = 16, ncol = 8)

gammas <- c(2, 5, 10 , 20)
horizons <-  c(1, 4, 8, 20, 40, 60, 80, 100)
k <- 1

for (i in c(1, 5, 9, 13)) {
  for (j in (1:length(horizons))) {
    BIG_TABLE[i:(i+3), j] <- rowMeans(port_calc(gamma = gammas[k], K = horizons[j])$alphafixedT)
    BIG_R_TABLE[i:(i+3), j] <- rowMeans(port_calc(gamma = gammas[k], K = horizons[j])$res_alphafixedT)
  }
  
  k <- k+1
}

HOR_ANA_VAL <- c(1, 3:100)

HOR_ANA <- matrix(0, nrow = 4, ncol = length(HOR_ANA_VAL))

for (i in (1:length(HOR_ANA_VAL))) {
  HOR_ANA[, i] <- rowMeans(port_calc(gamma = 5, K = HOR_ANA_VAL[i])$alphafixedT)
}

HOR_ANA_VAL <- c(1, 3:100)

HOR_ANA_R <- matrix(0, nrow = 4, ncol = length(HOR_ANA_VAL))

for (i in (1:length(HOR_ANA_VAL))) {
  HOR_ANA_R[, i] <- rowMeans(port_calc(gamma = 5, K = HOR_ANA_VAL[i])$res_alphafixedT)
}

GAM_ANA_VAL <- seq(1, 50, 0.05)

GAM_ANA <- matrix(0, nrow = 4, ncol = length(GAM_ANA_VAL))

for (i in (1:length(GAM_ANA_VAL))) {
  GAM_ANA[, i] <- rowMeans(port_calc(gamma = GAM_ANA_VAL[i], K = 100)$alphafixedT)
}

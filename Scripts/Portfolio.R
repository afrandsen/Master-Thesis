coef

Phi0_temp <- coef[ , ncol(coef)]
Phi1_temp <- coef[ , 1:(ncol(coef)-1)]

gamma <- 5
K <- 100
n <- 4

z <- data_sub[ , 1:(ncol(data_sub))]
z <- as.matrix(z)

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
  print('A0:')
  print(A0)
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


z2 <- t(z)
alphafixedT <- matrix(0, nrow = n-1, ncol = ncol(z2))

for (t in (2:ncol(z2))) {
  alphafixedT[, t] <- A0 + A1 %*% z2[, t-1]
  print(alphafixedT[, t])
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

alphafixedT <- rbind(alphafixedT, 1-colSums(alphafixedT))

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
  print(res_alphafixedT[, t])
}

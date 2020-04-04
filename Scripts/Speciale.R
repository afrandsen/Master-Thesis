library(IntroCompFinR)

set.seed(1)

mu_a <- 0.0225
mu_b <- 0.0394
mu_c <- 0.0128

rf <- 0.00125

return <- rbind(mu_a, mu_b, mu_c)

std_a <- 0.10347
std_b <- 0.12899
std_c <- 0.08660

std <- rbind(std_a, std_b, std_c)

corr <- matrix(c(1, 0.500, 0.500, 0.500, 1, 0.500, 0.500, 0.500, 1), nrow = 3)
cov <- sweep(sweep(corr, 1L, c(std_a,std_b,std_c), "*"), 2L, c(std_a,std_b,std_c), "*")


calc.portfolio <- function(r, cov, rf, a = c(1, 4, 8)){
  port <- list()
  
  gmv <- globalMin.portfolio(r, cov)
  
  tan <- tangency.portfolio(r, cov, rf)
  
  port[['gmv']] <- gmv
  port[['tan']] <- tan
  
  cara    <- matrix(0, nrow = length(a), ncol = length(c(r, rf)))
  cara_r  <- c()
  cara_sd <- c()
  
  for (i in seq_along(a)) {
    cara[i, ]  <- c(tan$weights*(tan$er-rf)/(a[i]*tan$sd^2), 1-(tan$er-rf)/(a[i]*tan$sd^2))
    
    cara_r[i]  <- sum(cara[i, ] * c(r, rf))
    
    cara_sd[i] <- sqrt(t(cara[i, 1:length(r)]) %*% cov %*% cara[i, 1:length(r)])
    
    port[[paste('cara', a[i], sep = '_')]] <- list("er" = as.vector(cara_r[i]),
                                                   "sd" = as.vector(cara_sd[i]),
                                                   "weights" = cara[i,])
  }
  
  return(port)
  
}


sim <- function(mu_a, mu_b, mu_c, std_a, std_b, std_c){
  e1 <- rnorm(1)
  e2 <- rnorm(1)
  e3 <- rnorm(1)
  
  r_a <- mu_a+std_a*e1
  r_b <- mu_b+std_b*(corr[, 2][1]*e1+sqrt(1-corr[, 2][1]^2)*e2)
  r_c <- mu_c + std_c * (corr[, 1][3] * e1 + ( (corr[, 2][3] - corr[, 1][2] * corr[, 1][3]) / (sqrt(1 - corr[, 1][2]^2)) ) * e2 + sqrt(1 - corr[, 1][3]^2 - ( (corr[, 2][3] - corr[, 1][2] * corr[, 1][3])^2 / (1- corr[, 1][2]^2) ) ) * e3)
  
  return(c(r_a, r_b, r_c))
}

data_1 <- t(matrix(replicate(200, sim(mu_a, mu_b, mu_c, std_a, std_b, std_c)), nrow = 3))
data_2 <- t(matrix(replicate(200, sim(mu_a, mu_b, mu_c, std_a, std_b, std_c)), nrow = 3))
data_3 <- t(matrix(replicate(200, sim(mu_a, mu_b, mu_c, std_a, std_b, std_c)), nrow = 3))

data_1_mean <- colMeans(data_1)
data_1_sd <- apply(data_1, 2, sd)
data_1_cor <- cor(data_1)
data_1_cov <- cov(data_1)

data_2_mean <- colMeans(data_2)
data_2_sd <- apply(data_2, 2, sd)
data_2_cor <- cor(data_2)
data_2_cov <- cov(data_2)

data_3_mean <- colMeans(data_3)
data_3_sd <- apply(data_3, 2, sd)
data_3_cor <- cor(data_3)
data_3_cov <- cov(data_3)

port_true <- calc.portfolio(return, cov, rf)
port_1    <- calc.portfolio(data_1_mean, data_1_cov, rf)
port_2    <- calc.portfolio(data_2_mean, data_2_cov, rf)
port_3    <- calc.portfolio(data_3_mean, data_3_cov, rf)





combine <- function(port){
  r <- c(port$gmv$er, port$tan$er, port$cara_1$er, port$cara_4$er, port$cara_8$er)
  sd <- c(port$gmv$sd, port$tan$sd, port$cara_1$sd, port$cara_4$sd, port$cara_8$sd)
  
  return(matrix(c(r, sd), ncol = 2))
}

r <- combine(port_true)[, 1]
r_1 <- combine(port_1)[, 1]
r_2 <- combine(port_2)[, 1]
r_3 <- combine(port_3)[, 1]

sd <- combine(port_true)[, 2]
sd_1 <- combine(port_1)[, 2]
sd_2 <- combine(port_2)[, 2]
sd_3 <- combine(port_3)[, 2]

# par(mfrow=c(1,1))
# plot(sort(sd), sort(r), type ='l', xlim=c(0,max(c(sd, sd_1, sd_2, sd_3))), ylim=c(0,max(c(r, r_1, r_2, r_3))))
# points(sort(sd), sort(r))
# lines(sort(sd_3), sort(r_3), col='red')
# points(sort(sd_3), sort(r_3), col='red')
# lines(sort(sd_2), sort(r_2), col='magenta')
# points(sort(sd_2), sort(r_2), col='magenta')
# lines(sort(sd_1), sort(r_1), col='blue')
# points(sort(sd_1), sort(r_1), col='blue')
# 
# par(mfrow=c(4,5))
# barplot(port_true$gmv$weights, xlab = 'Assets', ylab = 'Weight', main = 'GMV of TRUE')
# barplot(port_true$tan$weights, xlab = 'Assets', ylab = 'Weight', main = 'TAN of TRUE')
# barplot(port_true$cara_1$weights, xlab = 'Assets', ylab = 'Weight', main = 'a=1 of TRUE')
# barplot(port_true$cara_4$weights, xlab = 'Assets', ylab = 'Weight', main = 'a=4 of TRUE')
# barplot(port_true$cara_8$weights, xlab = 'Assets', ylab = 'Weight', main = 'a=8 of TRUE')
# 
# barplot(port_1$gmv$weights, xlab = 'Assets', ylab = 'Weight', main = 'GMV of 1')
# barplot(port_1$tan$weights, xlab = 'Assets', ylab = 'Weight', main = 'TAN of 1')
# barplot(port_1$cara_1$weights, xlab = 'Assets', ylab = 'Weight', main = 'a=1 of 1')
# barplot(port_1$cara_4$weights, xlab = 'Assets', ylab = 'Weight', main = 'a=4 of 1')
# barplot(port_1$cara_8$weights, xlab = 'Assets', ylab = 'Weight', main = 'a=8 of 1')
# 
# barplot(port_2$gmv$weights, xlab = 'Assets', ylab = 'Weight', main = 'GMV of 2')
# barplot(port_2$tan$weights, xlab = 'Assets', ylab = 'Weight', main = 'TAN of 2')
# barplot(port_2$cara_1$weights, xlab = 'Assets', ylab = 'Weight', main = 'a=1 of 2')
# barplot(port_2$cara_4$weights, xlab = 'Assets', ylab = 'Weight', main = 'a=4 of 2')
# barplot(port_2$cara_8$weights, xlab = 'Assets', ylab = 'Weight', main = 'a=8 of 2')
# 
# barplot(port_3$gmv$weights, xlab = 'Assets', ylab = 'Weight', main = 'GMV of 3')
# barplot(port_3$tan$weights, xlab = 'Assets', ylab = 'Weight', main = 'TAN of 3')
# barplot(port_3$cara_1$weights, xlab = 'Assets', ylab = 'Weight', main = 'a=1 of 3')
# barplot(port_3$cara_4$weights, xlab = 'Assets', ylab = 'Weight', main = 'a=4 of 3')
# barplot(port_3$cara_8$weights, xlab = 'Assets', ylab = 'Weight', main = 'a=8 of 3')



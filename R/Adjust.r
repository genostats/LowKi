v <- function(P1, P2) P1*(1-P1) + 4*P2*(1-P2) - 4*P1*P2
Adjust <- function(set_s, X, gl.Aa, gl.aa) {
  ## for each pair selected, run our linear regression
  f_t <- NULL
  for(m in 1:nrow(set_s)){
    i <- set_s[m, 1];
    j <- set_s[m, 2]
    
    ## calculate point wise estimates of kinship phiA
    phi <- X[i,]*X[j,]
    v1 <- v(gl.Aa[i,], gl.aa[i,])
    v2 <- v(gl.Aa[j,], gl.aa[j,])
    ## Calculate the mean of the 'fuzziness' of the two genotypes
    ## depending on whether we are interrogating a diagonal or off-diagonal
    ## pair, run one of two linear regression models and keep the intercept
    ## estimate.
    interc <- lm(phi ~ v1*v2)$coefficients[1]
    f_t <- rbind(f_t, c(interc, mean(phi, na.rm=TRUE)))
  }
  
  ### We need to compare the estimates of off-diagonal elements to the raw
  ### estimates, and the diagonal elements to the assumed value of 1.00
  w <- which(set_s[,1] == set_s[,2])
  B1 <- lm(f_t[-w,1] ~ f_t[-w,2])$coefficients[1:2]
  if( abs(B1[1]) > 0.005 ){ f_t[,1] <- f_t[,1] - B1[1] }
  beta1 <- B1[2]
  beta2 <- mean(1/f_t[w,1])
  ## Combining this information gives a reasonable beta
  beta <- beta1*beta2
  list(beta = beta1*beta2, beta1 = beta1)
}

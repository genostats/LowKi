#' Title
#'
#' @param Glmat 
#' @param method 
#' @param userFun 
#' @param adjust 
#' @param adjust.param 
#'
#' @return
#' @export
#'
#' @examples
lowGRM <- function(Glmat, method = "phred", userFun, adjust = TRUE, adjust.param = c(25,50)) {
  
  R <- rescale.GL(Glmat, method, userFun)
  gL1AA <- R[[1]]
  gL1Aa <- R[[2]]
  gL1aa <- R[[3]]
  
  ## Reweight to give three genotype probabilities
  
  ## Estimate the population genotype frequencies
  LAA <- colMeans(gL1AA, na.rm=T)
  LAa <- colMeans(gL1Aa, na.rm=T)
  Laa <- colMeans(gL1aa, na.rm=T)
  
  ## Estimate the varaince of genotypes in the population
  Vad <- LAa + (4*LAA*Laa) - LAa*LAa
  
  ## Additive components described in the paper
  UA_1 <- (1/sqrt(Vad))*(-LAa-2*Laa);
  UA_2 <- (1/sqrt(Vad))*(1-LAa-2*Laa);
  UA_3 <- (1/sqrt(Vad))*(2-LAa-2*Laa);

  UA_1[which(is.infinite(UA_1))] <- NA
  UA_2[which(is.infinite(UA_2))] <- NA
  UA_3[which(is.infinite(UA_3))] <- NA
  
  
  ## Sweep these across across the genotype probability matrices and sum up
  XA   <- sweep(gL1AA, 2, UA_1, "*")
  XA_2 <- sweep(gL1Aa, 2, UA_2, "*")
  XA_3 <- sweep(gL1aa, 2, UA_3, "*")
  XA   <- as.matrix(XA + XA_2 + XA_3) 
  

  ## Compute the matrix
  cat("computing matrix...\n")
  
  K_GL_A <- mmult(XA)
  ## Give the matrix col and rownames of individuals.
  rownames(K_GL_A) <- Glmat[[4]];
  colnames(K_GL_A) <- Glmat[[4]]
  
  ## Readjustment method...
  if(adjust) {
    cat("adjusting...\n")
    
    ## sample 25 diagonal elements at random
    s1 <- sample(1:nrow(gL1AA), adjust.param[1], replace = FALSE)
    set_s <- cbind(s1,s1)
    colnames(set_s) <- c("i","j")
    A1 <- reshape.GRM(K_GL_A);
    A1 <- arrange(A1, by=k)
    ## select the 10 smallest, 10 biggest, and 30 random off-diagonal values (or user defined values).
    s2 <- c(1:10, sample(11:(nrow(A1)-10), max(adjust.param[2]-20,10), replace = FALSE), (nrow(A1)-9):nrow(A1))
    
    set_s <- rbind(set_s, A1[s2,1:2])
  
    # set_s <- rbind( cbind(1:nrow(XA), 1:nrow(XA)), cbind(A1$i, A1$j) )
    b <- Adjust(set_s, XA, gL1Aa, gL1aa)
   
    # adjust diagonals and off-diagonals seperately.
    K_GL_A <- b$beta*K_GL_A;
    diag(K_GL_A) <- diag(K_GL_A)/b$beta1
  }
  return(K_GL_A)
}


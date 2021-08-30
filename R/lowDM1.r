lowDM1 <- function(a, adjust = TRUE, adjust.param = c(25,50)) {
  
  gL1AA <- 1-a$P1-a$P2
  gL1Aa <- a$P1
  gL1aa <- a$P2

  LAA<-colMeans(gL1AA,na.rm=T)
  LAa<-colMeans(gL1Aa,na.rm=T)
  Laa<-colMeans(gL1aa,na.rm=T)
  
  Gamma <- Laa + ((4*Laa*LAA)/LAa) + LAA
  Gamma <- 1/(sqrt(Gamma))
  
  UD_1 <- Gamma*sqrt(Laa/LAA);
  UD_2 <- (-2)*Gamma*sqrt((Laa*LAA)/(LAa*LAa));
  UD_3 <- Gamma*sqrt(LAA/Laa);

  UD_1[which(is.infinite(UD_1))] <- NA
  UD_2[which(is.infinite(UD_2))] <- NA
  UD_3[which(is.infinite(UD_3))] < -NA
  
  XD   <- sweep(gL1AA, 2, UD_1, "*")
  XD_2 <- sweep(gL1Aa, 2, UD_2, "*")
  XD_3 <- sweep(gL1aa, 2, UD_3, "*")
  
  XD <- as.matrix(XD + XD_2 + XD_3)

  cat("computing matrix...\n")

  D_GL_D <- mmult(XD)
  rownames(D_GL_D) <- Glmat[[4]];
  colnames(D_GL_D) <- Glmat[[4]]
  
  if(adjust) {
    cat("adjusting...\n")
   
 
    s1 <- sample(1:nrow(gL1AA), adjust.param[1], replace = FALSE)
    set_s <- cbind(s1,s1)
    colnames(set_s) <- c("i","j")
    A1 <- reshape.GRM(D_GL_D);
    A1 <- arrange(A1, by=k)
    
    s2 <- c(1:10, sample(11:(nrow(A1)-10), max(adjust.param[2]-20,10), replace = FALSE), (nrow(A1)-9):nrow(A1))
    
    set_s <- rbind(set_s, A1[s2,1:2])
   
    b <- Adjust(set_s, XD, gL1Aa, gL1aa)


    D_GL_D <- b$beta*D_GL_D;
    diag(D_GL_D) <- diag(D_GL_D)/b$beta1
  }
  return(D_GL_D)
}

rescale.GL <- function(Glmat, method = "none" , userFun) {
  ## Create NxM matricies for the three genotype likelihoods, rescaled from an
  ## anticipated phred likelihood calculation from GATK.
  if (method=="none") {
    gL1AA <- t(Glmat[[1]])
    gL1Aa <- t(Glmat[[2]]) 
    gL1aa <- t(Glmat[[3]]) 
  } else if(method == "phred"){
    ## Assuming its a VCF with phred likelihoods -  PL = -10*log10(X) - K for
    ## some constant K which we can ignore as on the log scale and we will
    ## re-weight later.
    ## change
    f <- function(x) 10^(-0.1*as.numeric(as.character(x)))
    gL1AA <- t(mutate_all(Glmat[[1]], f))
    gL1Aa <- t(mutate_all(Glmat[[2]], f)) 
    gL1aa <- t(mutate_all(Glmat[[3]], f))
  } else if(method == "imputation"){
    f <- function(x) as.numeric(as.character(x))
    gL1AA <- t(mutate_all(Glmat[[1]], f))
    gL1Aa <- t(mutate_all(Glmat[[2]], f))
    gL1aa <- t(mutate_all(Glmat[[3]], f))
  } else if (method=="custom"){
    ## Perhaps the user has a differnt type of likelihood (from an alternate
    ## software and so could specify their own function. In particular if
    ## Genotype Likelihoods have already been calculated as in BEAGLE
    gL1AA <- t(mutate_all(Glmat[[1]], userFun))
    gL1Aa <- t(mutate_all(Glmat[[2]], userFun))
    gL1aa <- t(mutate_all(Glmat[[3]], userFun))
  } else {
    stop("method ???")
  }
  weight <- gL1AA + gL1Aa + gL1aa;
  gL1AA  <- gL1AA / weight;
  gL1Aa  <- gL1Aa / weight;
  gL1aa  <- gL1aa / weight
  list( gL1AA, gL1Aa, gL1aa )
}

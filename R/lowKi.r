#' Computing Kinship and Fraternity matrices
#'
#' @param filename path to a VCF file
#' @param field which field of the VCF file to use 
#' @param freqs [optional] allele frequencies (see details)
#' @param adjust logical. TRUE to use adjustment method
#' @param fraternity logical. TRUE to compute the fraternity matrix instead of the kinship matrix
#' @param adjust.par Parameters for the adjustment procedure (details below).
#'
#' @details This function computes the kinship matrix or the fraternity matrix (aka as 'dominance matrix')
#' from shallow sequencing data. It relies on moment estimates and on an adjustment method which allows
#' to correct the bias.
#' @details The method is exposed in details in \url{https://doi.org/10.21203/rs.3.rs-1109592/v1}. 
#' The three values in `adjust.par`
#' allow to control the number of individuals used in the adjustment. The first (resp. second) value gives the number 
#' of pairs of individuals with the lowest (resp. highest) kinship to be selected; the third gives the number of
#' random individuals to be added. The adjustment procedure is based on all the estimated coefficients between 
#' the selected individuals. Note that the two first values should be at least 5; in case they are lower, 
#' they are replaced by 5. Using `adjust.par = c(0,0,n)` with `n` being the size of the sample results in using
#' all pairs of individuals in the adjustment procedure, which is computationnally more heavy.
#' @details The argument `freqs` is optionnal. If provided, it should be a data frame which a column `p` giving
#' the alternate allele frequency of each SNP in the vcf file.
#' @details All SNPs present in the VCF file will be used, regardless of the MAF or of the chromosome number.
#'
#' @return A symmetric matrix.
#' @export
#'
#' @examples # see vignette for more examples
#' vcf.file <- system.file("extdata", "shallow.vcf.gz", package="LowKi")
#' kinship.file <-system.file("extdata", "kinship.rds", package="LowKi")
#' K.low.1 <- lowKi(vcf.file)
#' # Another solution: using MLE allele frequencies
#' MLfreq <- vcf.allele.freq(vcf.file)
#' K.low.2 <- lowKi(vcf.file, freqs = MLfreq)
#' K.real <- readRDS(kinship.file)
#' par(mfrow=c(1,2))
#' plot(K.real, K.low.1, xlab = "true values", ylab = "LowKi estimates")
#' abline(0,1,col="red")
#' plot(K.real, K.low.2, xlab = "true values", ylab = "LowKi estimates")
#' abline(0,1,col="red")



lowKi <- function(filename, field = c("PL", "GP"), freqs, adjust = TRUE, fraternity = FALSE, adjust.par = c(20,20,100) ) {

  filename <- path.expand(filename)
  field <- match.arg(field)

  # matrix of unadjusted coefficients
  if(missing(freqs)) {
    K <- RawKinVcfLoki(filename, field, fraternity)
  } else {
    K <- RawKinVcfClassic(filename, field, freqs$p, fraternity)
  }

  if(adjust) {
    # extracting indices of extreme off diagonal values from the matrix
    # this could be rewritten in C++ to speed up a little bit
    # but I don't think it's a bottleneck
    n <- nrow(K)
    o <- order(K) - 1
    o <- cbind(o %% n, o %/% n) + 1
    o <- o[ o[,1] != o[,2] , ]
  
    n.lo <- max(5, adjust.par[1])
    n.hi <- max(5, adjust.par[2])
    n.random <- adjust.par[3]
    s <- c(as.vector(head(o, n.lo)), as.vector(tail(o, n.hi)))
  
    # plus a few random indices
    if(n.random < n) {
      s <- c(s, sample(1:n, n.random))
      s <- sort(unique(s))
    } else {
      s <- 1:n
    }

    # computing regression adjusted matrix for all pairs of individuals in s
    if(missing(freqs)) {
      L <- PartialKinVcfLoki(filename, s - 1L, field, TRUE, fraternity, TRUE) 
    } else {
      L <- PartialKinVcfClassic(filename, s - 1L, field, freqs$p, TRUE, fraternity, TRUE) 
    }
    # estimating beta1 and beta2 multiplicative coeffs
    off.diag.raw <- K[s,s][upper.tri(L)]
    off.diag.adj <- L[ upper.tri(L) ]
    B1 <- lm(off.diag.adj ~ off.diag.raw)$coefficients[1:2]

    # adjusting K
    alpha <- B1[1] 
    beta1 <- B1[2]
    beta2 <- mean( 1/(diag(L) - alpha) )

    diag.K <- diag(K)
    K <- beta1 * beta2 * K

    # final shift using 25th percentile
    diag(K) <- NA
    K <- K - quantile(K, 0.25, na.rm=TRUE)

    # adjusting the diagonal
    if(fraternity) {
      diag(K) <- 1
    } else { 
      diag(K) <- diag.K * beta2
    }
  }
  K
}

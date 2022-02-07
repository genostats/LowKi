#' Computing Kinship and Fraternity matrices
#'
#' @param filename path to a VCF file
#' @param field which field of the VCF file to use 
#' @param adjust logical. TRUE to use adjustment method
#' @param fraternity logical. TRUE to compute the fraternity matrix instead of the kinship matrix
#' @param adjust.par Parameters for the adjustment procedure (details below).
#'
#' @details This function computes the kinship matrix or the fraternity matrix (aka as 'dominance matrix')
#' from shallow sequencing data. It relies on moment estimates and on an adjustment method which allows
#' to correct the bias.
#' @details The method is exposed in details in a yet unpublished paper. The three values in `adjust.par`
#' allow to control the number of individuals used in the adjustment. The first (resp. second) value gives the number 
#' of pairs of individuals with the lowest (resp. highest) kinship to be selected; the third gives the number of
#' random individuals to be added. The adjustment procedure is based on all the estimated coefficients between 
#' the selected individuals.
#' @details All SNPs present in the VCF file will be used, regardless of the MAF or of the chromosome number.
#'
#' @return A symmetric matrix.
#' @export
#'
#' @examples # see vignette for more examples
#' vcf.file <- system.file("extdata", "shallow.vcf.gz", package="LowKi")
#' kinship.file <-system.file("extdata", "kinship.rds", package="LowKi")
#' K.low <- lowKi(vcf.file)
#' K.real <- readRDS(kinship.file)
#' plot(K.real, K.low, xlab = "true values", ylab = "LowKi estimates")
#' abline(0,1,col="red")



lowKi <- function(filename, field = c("PL", "GP"), adjust = TRUE, fraternity = FALSE, adjust.par = c(20,20,10) ) {

  filename <- path.expand(filename)
  field <- match.arg(field)

  # matrix of unadjusted coefficients
  K <- RawKinVcf(filename, field, fraternity)
 
  if(adjust) {
    # extracting indices of extreme off diagonal values from the matrix
    # this could be rewritten in C++ to speed up a little bit
    # but I don't think it's a bottleneck
    n <- nrow(K)
    o <- order(K) - 1
    o <- cbind(o %% n, o %/% n) + 1
    o <- o[ o[,1] != o[,2] , ]
  
    n.lo <- adjust.par[1]
    n.hi <- adjust.par[2]
    n.random <- adjust.par[3]
    s <- c(as.vector(head(o, n.lo)), as.vector(tail(o, n.hi)))
  
    # plus a few random indices
    s <- c(s, sample(1:n, n.random))
    s <- sort(unique(s))
  
    # computing regression adjusted matrix for all pairs of individuals in s
    L <- PartialKinVcf(filename, s - 1L, field, TRUE, fraternity, TRUE) 
  
    # estimating beta1 and beta2 multiplicative coeffs
    off.diag.raw <- K[s,s][upper.tri(L)]
    off.diag.adj <- L[ upper.tri(L) ]
    B1 <- lm(off.diag.adj ~ off.diag.raw)$coefficients[1:2]
  
    beta1 <- B1[2]
    beta2 <- mean( 1/(diag(L) - B1[1]) )
  
    # adjusting K
    diag.K <- diag(K)
    K <- K * beta1 * beta2
    diag(K) <- diag.K * beta2
  }

  K
}

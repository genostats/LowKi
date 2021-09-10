#' Title
#'
#' @param filename path to a VCF file
#' @param field which field of the VCF file to use 
#' @param adjust
#' @param dominance
#' @param constraint constraint regression parameters of v1 and v2 to be equal
#'
#' @details I left the 'constraint' option for the moment but this should disappear
#' @return
#' @export
#'
#' @examples

lowKin.vcf <- function(filename, field = c("PL", "GP"), adjust = TRUE, dominance = FALSE, adjust.par = c(20,20,10) ) {

  filename <- path.expand(filename)
  field <- match.arg(field)

  # matrix of unadjusted coefficients
  K <- RawKinVcf(filename, field, dominance)

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
  s <- unique(s)

  # computing regression adjusted matrix for all pairs of individuals in s
  L <- PartialKinVcf(filename, s - 1L, field, TRUE, dominance, TRUE) 

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
  K
}

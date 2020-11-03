#' Title
#'
#' @param Glmat 
#' @param method 
#' @param adjust 
#' @param constraint constraint regression parameters of v1 and v2 to be equal
#'
#' @return
#' @export
#'
#' @examples
lowKin <- function(Glmat, method = c("none","imputation", "phred"), adjust = TRUE, dominance = FALSE, constraint = FALSE) {
  method <- match.arg(method)
  R <- rescale.GL(Glmat, method)
  gL1Aa <- R[[2]]
  gL1aa <- R[[3]]
  K <- lowKincpp(R[[2]], R[[3]], adjust, dominance, constraint)
  # we could do that even when adjust is FALSE...
  if(adjust) 
    K/mean(diag(K))
  else
    K
}

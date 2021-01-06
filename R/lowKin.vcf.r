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

lowKin.vcf <- function(filename, field = c("PL", "GP"), adjust = TRUE, dominance = FALSE, constraint = TRUE) {
  filename <- path.expand(filename)
  field <- match.arg(field)
  lowKinVcf(filename, field, adjust, dominance, constraint)
}


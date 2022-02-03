#' Maximum likelihood estimates of allele frequencies
#'
#' @param filename path to a VCF file
#' @details Estimates allele frequencies by maximizing a likelihood on data from 'AD' field.
#' @return a data frame with columns `snp.id`, `chr`, `pos`, `A1`, `A2`, `p`. The last column
#' is the frequency of the alternate allele `A2`.
#' @export
#'
#' @examples # computing allele freqs in example file
#' vcf.file <- system.file("extdata", "shallow.vcf.gz", package="LowKi")
#' MLfreq <- vcf.allele.freq(vcf.file)
#' # comparison of allele freqs from hardcalled genotypes
#' z <- read.vcf(vcf.file)
#' plot(MLfreq$p, z@p, pch = ".", cex = 2)


vcf.allele.freq <- function(filename, field = c("AD", "PL", "GP")) {
  field <- match.arg(field) 
  if(field == "AD") 
    data.frame(vcfAlleleFreqAD(filename))
  else 
    data.frame(vcfAlleleFreqPr(filename, field))
}

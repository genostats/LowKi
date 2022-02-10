require(LowKi)
vcf.file <- system.file("extdata", "shallow.vcf.gz", package="LowKi")
vcf.file.20 <- system.file("extdata", "shallow_20.vcf.gz", package="LowKi")
kinship.file <-system.file("extdata", "kinship.rds", package="LowKi")
kinship.file.20 <-system.file("extdata", "kinship_20.rds", package="LowKi")

MLfreq <- vcf.allele.freq(vcf.file)

K.low <- lowKi(vcf.file, freqs = MLfreq, adjust = FALSE)
K.low.20 <- lowKi(vcf.file.20, freqs = MLfreq, adjust = FALSE)

K.real <- readRDS(kinship.file)
K.real.20 <- readRDS(kinship.file.20)

par(ask = TRUE)
m <- match( rownames(K.low.20), rownames(K.low) )
# raw coeffs with same allelic ref are identical, as they ought to be
plot( K.low.20, K.low[m,m]) ; abline(0,1,col="red")

K.low.adj <- lowKi(vcf.file, freqs = MLfreq, adjust = TRUE)
K.low.adj.20 <- lowKi(vcf.file.20, freqs = MLfreq, adjust = TRUE)
plot( K.low.adj.20, K.low.adj[m,m]); abline(0,1,col="red")

K.low.adj.noref <- lowKi(vcf.file, adjust = TRUE)
K.low.adj.noref.20 <- lowKi(vcf.file.20, adjust = TRUE)
plot( K.low.adj.noref.20, K.low.adj.noref[m,m]); abline(0,1,col="red")


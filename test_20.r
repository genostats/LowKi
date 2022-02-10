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

m <- match( rownames(K.low.20), rownames(K.low) )
plot( K.low.20, K.low[m,m])

K.low.adj <- lowKi(vcf.file, freqs = MLfreq, adjust = TRUE)
K.low.adj.20 <- lowKi(vcf.file.20, freqs = MLfreq, adjust = TRUE)

plot( K.low.adj.20, K.low.adj[m,m])


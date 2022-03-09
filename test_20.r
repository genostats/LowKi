require(LowKi)
vcf.file <- system.file("extdata", "shallow.vcf.gz", package="LowKi")
vcf.file.20 <- system.file("extdata", "shallow_20.vcf.gz", package="LowKi")
kinship.file <-system.file("extdata", "kinship.rds", package="LowKi")
kinship.file.20 <-system.file("extdata", "kinship_20.rds", package="LowKi")

fraternity.file.20 <- system.file("extdata", "fraternity_20.rds", package="LowKi")
af.file <- system.file("extdata", "alleleFreqs.tsv", package="LowKi")

# MLfreq <- vcf.allele.freq(vcf.file)
AF <- read.table(af.file, header = TRUE)

# K.low <- lowKi(vcf.file, freqs = MLfreq, adjust = FALSE)
K.low.20 <- lowKi(vcf.file.20)
K.low.20.af <- lowKi(vcf.file.20, freqs = AF)

D.low.20 <- lowKi(vcf.file.20, fraternity = TRUE)
D.low.20.af <- lowKi(vcf.file.20, freqs = AF, fraternity = TRUE)


# K.real <- readRDS(kinship.file)
K.real.20 <- readRDS(kinship.file.20)
D.real.20 <- readRDS(fraternity.file.20)

par(mfrow=c(2,2))
plot(K.real.20, K.low.20); abline(0,1,col="red")
plot(D.real.20, D.low.20); abline(0,1,col="red")

plot(K.real.20, K.low.20.af); abline(0,1,col="red")
plot(D.real.20, D.low.20.af); abline(0,1,col="red")



stop()
m <- match( rownames(K.low.20), rownames(K.low) )
# raw coeffs with same allelic ref are identical, as they ought to be
plot( K.low.20, K.low[m,m]) ; abline(0,1,col="red")

K.low.adj <- lowKi(vcf.file, freqs = MLfreq, adjust = TRUE)
K.low.adj.20 <- lowKi(vcf.file.20, freqs = MLfreq, adjust = TRUE)
plot( K.low.adj.20, K.low.adj[m,m]); abline(0,1,col="red")

K.low.adj.noref <- lowKi(vcf.file, adjust = TRUE)
K.low.adj.noref.20 <- lowKi(vcf.file.20, adjust = TRUE)
plot( K.low.adj.noref.20, K.low.adj.noref[m,m]); abline(0,1,col="red")


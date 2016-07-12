(N <- 500)
(MU <- c(0,0))
(SIGMA <- matrix(c(1,0,0,1), nrow = 2))



sample = rnormal(n = N, mu = c(1,2), sigma = SIGMA)
sample

apply(sample, 2, mean)
cov(sample)

set.seed(1)
sample = rnormalmultinomial(size= rep(5, N), mu = MU, sigma = SIGMA)
sample


apply(sample, 2, mean)
cov(sample)

df = as.data.frame(sample)
barplot(apply(df, 2, sum))

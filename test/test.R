(N <- 500)
(MU <- c(0,0))

SIGMA.ilr <- matrix(c(1,0,
                       0,1), nrow = 2)
ilr_to_alr = MASS::ginv(ilrBase(D = 3)) %*% matrix(c(1,0,-1,
                                                     0,1,-1), ncol=2)
SIGMA = t(ilr_to_alr) %*% SIGMA.ilr %*% ilr_to_alr

set.seed(1)
sample_normal = rnormal(n = N, mu = c(1,2), sigma = SIGMA)

set.seed(1)
sample_multinomial = rmultinomial(n = N, size = c(10, 100), prob = c(0.2,0.8))

set.seed(1)
sample_nm = rnormalmultinomial(size= rep(5, N), mu = MU, sigma = SIGMA, probs = TRUE)

library(ggtern)
probs = as.data.frame(sample_nm$probs)
ggtern() +
  geom_point(data=probs, aes(x = V1, y = V2, z = V3))

plot(as.data.frame(ilr(probs)), asp = 1)

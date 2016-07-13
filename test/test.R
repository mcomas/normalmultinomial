library(mixpack)

(N <- 10000)
(MU <- c(0,0))

SIGMA.ilr <- matrix(c(1,0,
                       0,1), nrow = 2)
ilrB = t(clr_coordinates(t(matrix(unlist(ilr_basis(3)), ncol=2))))
ilr_to_alr = MASS::ginv(ilrB) %*% matrix(c(1,0,-1,
                                                     0,1,-1), ncol=2)
SIGMA = t(ilr_to_alr) %*% SIGMA.ilr %*% ilr_to_alr

set.seed(1)
sample_normal = rnormal(n = N, mu = c(1,2), sigma = SIGMA)

set.seed(1)
sample_multinomial = rmultinomial(n = N, size = c(10, 100), prob = c(0.2,0.8))

set.seed(1)
sample_nm = rnormalmultinomial(N, 100, mu = MU, sigma = SIGMA, probs = TRUE)

library(ggtern)
probs = as.data.frame(sample_nm$probs)
counts = as.data.frame(sample_nm$counts)
ggtern() +
  geom_point(data=probs, aes(x = V1, y = V2, z = V3), col='red', alpha=0.05) +
  theme_bw()
ggtern() +
  geom_point(data=counts, aes(x = V1, y = V2, z = V3), col='red', alpha=0.015, size = 3) +
  theme_bw()

# plot(as.data.frame(ilr_coordinates(probs)), asp = 1)

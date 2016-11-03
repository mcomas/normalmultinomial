library(normalmultinomial)

p = c(0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
      0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
      0.050, 0.050, 0.050, 0.050, 0.500)
h = ilr_coordinates(p)
sigma = diag(length(h))

X = rnormalmultinomial(100, 20, h, sigma)
fit = nm_fit2(X, eps = 0.02)

p = c(0.25,0.25,0.5)
h = ilr_coordinates(p)
sigma = diag(length(h))

X = rnormalmultinomial(100, 20, h, sigma)
MU = ilr_coordinates(colSums(X/rowSums(X)))
SIGMA = sigma
Z = matrix(randtoolbox::halton(10000, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))

system.time(res1 <- apply(X, 1, expectedMonteCarlo, MU, SIGMA, Z))
system.time(res2 <- apply(X, 1, expectedMonteCarlo2, MU, SIGMA, Z))

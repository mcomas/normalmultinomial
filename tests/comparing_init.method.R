library(normalmultinomial)
library(parallel)

load('data/dataset-N_01000-n_00050-s_00001-seed_00010.RData')


fdm = dm_fit(XZ)
eps = 0.0005*length(p)

cl = makeCluster(3)
fnm1 <- nm_fit(XZ, verbose = TRUE, eps = eps, nsim = 100, parallel.cluster = cl, init.method = 'dm')
stopCluster(cl)

cl = makeCluster(3)
fnm2 <- nm_fit(XZ, verbose = TRUE, eps = eps, nsim = 100, parallel.cluster = cl, init.method = 'approx')
stopCluster(cl)

cl = makeCluster(3)
fnm3 <- nm_fit(XZ, verbose = TRUE, eps = eps, nsim = 100, parallel.cluster = cl, init.method = 'approx2')
stopCluster(cl)


compare = function(H1, H2) mean(apply(H1-H2, 1, function(x) sqrt(sum(x^2))))
p = params[[s]]$p
P.gs = matrix(p, nrow = N, ncol = length(p), byrow = TRUE)
H.gs = ilr_coordinates(P.gs)

P = fdm$expected
H = ilr_coordinates(P)
data.frame('paired.dist' = compare(H, H.gs))

P = fnm1$expected
H = ilr_coordinates(P)
data.frame('paired.dist' = compare(H, H.gs))

P = fnm2$expected
H = ilr_coordinates(P)
data.frame('paired.dist' = compare(H, H.gs))

P = fnm3$expected
H = ilr_coordinates(P)
data.frame('paired.dist' = compare(H, H.gs))




cl = makeCluster(3)
fit_nm2 <- nm_fit(XZ, verbose = TRUE, eps = eps, nsim = 100, parallel.cluster = cl)
stopCluster(cl)

cl = makeCluster(3)
fit_nm3 <- nm_fit(XZ, verbose = TRUE, eps = eps, nsim = 100, parallel.cluster = cl, init = 'sum')
stopCluster(cl)

library(parallel)
library(normalmultinomial)

if(!exists('build')) build = 'N_01000-n_00050-s_00001-seed_00001'
pattern_build = "N_([0-9]+)-n_([0-9]+)-s_([0-9]+)-seed_([0-9]+)*"

load(sprintf('datasets/dataset-%s.RData', build))

eps = 0.0005*length(p)

cl = makeCluster(3)
fit_nm1 <- nm_fit(XZ, verbose = TRUE, eps = eps, nsim = 100, parallel.cluster = cl, init = 'mean')
stopCluster(cl)

cl = makeCluster(3)
fit_nm2 <- nm_fit(XZ, verbose = TRUE, eps = eps, nsim = 100, parallel.cluster = cl)
stopCluster(cl)

cl = makeCluster(3)
fit_nm3 <- nm_fit(XZ, verbose = TRUE, eps = eps, nsim = 100, parallel.cluster = cl, init = 'sum')
stopCluster(cl)

save(fit, file = sprintf('datasets/replacement-%s-method_nm.RData', build))

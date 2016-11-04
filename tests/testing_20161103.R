library(normalmultinomial)

p = c(0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
      0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
      0.050, 0.050, 0.050, 0.050, 0.500)
h = ilr_coordinates(p)
sigma = diag(length(h))

X = rnormalmultinomial(100, 20, h, sigma)
fit = nm_fit3(X, eps = 0.02, verbose = TRUE)
fit = nm_fit2(X, eps = 0.02, verbose = TRUE)
fit = nm_fit(X, eps = 0.02, verbose = TRUE)

p = c(0.25,0.25,0.5)
h = ilr_coordinates(p)
sigma = diag(length(h))

X = rnormalmultinomial(100, 20, h, sigma)
MU = ilr_coordinates(colSums(X/rowSums(X)))
SIGMA = sigma
Z = matrix(randtoolbox::halton(100, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))

library(microbenchmark)
method_old = function() apply(X, 1, expectedMonteCarlo, MU, SIGMA, Z)
method_new = function() apply(X, 1, expectedMonteCarlo2, MU, SIGMA, Z)
method_new_E = function() sapply(1:nrow(X), function(i) expectedMonteCarlo3(X[i,], MU, SIGMA, Z, E[i,]), simplify = FALSE)

FIT = method_old()
E = t(sapply(FIT, function(fit) colMeans(na.omit(fit[[2]]))))

microbenchmark(
  method_old(),
  method_new(),
  method_new_E(), times = 500)

res1 = method_old()
res2 = method_new()
res3 = method_new_E()

head(res1[[1]][[2]])
head(res2[[1]][[2]])
head(res3[[1]][[2]])

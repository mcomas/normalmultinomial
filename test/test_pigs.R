library(mixpack)
load('data/Pigs.rda')

X = as.matrix(Pigs)

fit = normalmultinomial_fitting(X)
(MU <- fit[[1]])
(SIGMA <- fit[[2]])

A = stepE(X = X, mu = MU, sigma = SIGMA)

Prob = cbind(exp(A), 1) / apply(cbind(exp(A), 1), 1, sum)
X.replaced = apply(X, 1, sum) * Prob

biplot(princomp(clr_coordinates(X.replaced)), asp=1, scale=0, xlim = c(-5, 5), ylim = c(-5,5))

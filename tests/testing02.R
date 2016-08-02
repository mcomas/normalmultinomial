invalr = function(alr){
  a = exp(cbind(alr, 0))
  a/sum(a)
}
X = diag(3)
MU = matrix(c(0,0), 1)
SIGMA.ILR = diag(2)
SIGMA.ALR = matrix(c(2,1,1,2), nrow=2)

t(solve(ilr_to_alr(3)) %*% apply(X, 1, function(x) expectedA2(x, MU, SIGMA.ALR, nsim=10000)[[1]]))
expected_initial(matrix(X, ncol=3), MU, SIGMA.ILR, se_eps = 0.0005)[[1]][,2:3]



stepEM1(X, MU, SIGMA.ALR)
(it <- stepEM3(X, MU, SIGMA.ALR))

(it <- stepEM4(X, it[[1]], it[[2]], it[[3]], it[[4]], 10000))


solve(ilr_to_alr(3)) %*% expectedA2(X[1,], MU, SIGMA.ALR, nsim=10000)[[2]] %*% t(solve(ilr_to_alr(3)))
expected_initial(matrix(X[1,], ncol=3), MU, SIGMA.ILR, se_eps = 0.001)[[1]][,4:6]



fit = normalmultinomial_fitting(X, nsim = 100000)
a = exp(cbind(stepE(X, fit[[1]], fit[[2]], nsim = 100000),0))
a/apply(a,1,sum)

(E <- expectedA2(X, MU, SIGMA, se_eps = 0.001))
round(inv_ilr_coordinates(E[[1]][,2:3]),3)

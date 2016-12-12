library("MASS")
data(menarche)

initialize_with_dm = function(X){
  E = dm_fit(X)$expected

  H = ilr_coordinates(E)
  MU = colMeans(H)
  SIGMA = cov(H)
  list(H = H, MU = MU, SIGMA = SIGMA)
}
generate_mv_normal_rnd = function(n, dim){
  if(dim <= 6){
    Z = matrix(randtoolbox::halton(n, dim = dim, normal = TRUE), ncol=dim)
  }else{
    Z = matrix(randtoolbox::sobol(n, dim = dim, normal = TRUE), ncol=dim)
  }
  return(Z)
}

Y = with(menarche, cbind(Menarche, Total-Menarche))
X = cbind(1,menarche$Age)

init = initialize_with_dm(Y)
E = init$H
B = matrix(c(init$MU,
             rep(0, length(init$MU))), ncol = length(init$MU))

MU = X %*% B
SIGMA = init$SIGMA

Z = generate_mv_normal_rnd(10000, nrow(SIGMA))

pB = B
err = 1
eps = 10^-6
max.em.iter = 100
iter = 0
while(err > eps & iter < max.em.iter){
  fit = lapply(1:nrow(Y), function(i) expectedMonteCarloFast(Y[i,], MU[i,], SIGMA, Z, E[i,]))
  y = sapply(fit, function(f) f[[2]])
  m = lm(y~0+X)
  B = cbind(coef(m))
  err = max(B-pB)
  pB = B
  MU = cbind(m$fitted.values)
  SIGMA = matrix(var(m$residuals))
  iter = iter + 1
}

summary(m)

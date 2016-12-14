library(normalmultinomial)

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

Y = as.matrix(parliament2015[,c('jxsi', 'psc', 'pp', 'catsp', 'cs', 'cup')])
X = cbind(1,
          boot::logit(parliament2015$birth.cat/parliament2015$pop),
          rnorm(nrow(parliament2015)) )
# parliament2015$ss.agri,
# parliament2015$lat,
# parliament2015$long
init = initialize_with_dm(Y)
E = init$H
B = t(apply(E, 2, function(y){
  b = solve(t(X) %*% X) %*% (t(X) %*% y)
}))

MU = X %*% t(B)
SIGMA = diag(mean(diag(init$SIGMA)), ncol(E))

Z = generate_mv_normal_rnd(1000, nrow(SIGMA))

pB = B
err = 1
eps = 0.001
max.em.iter = 100
iter = 0
while(err > eps & iter < max.em.iter){

  FIT = lapply(1:nrow(Y), function(i) expectedMonteCarloFast(Y[i,], MU[i,], SIGMA, Z, E[i,]))
  E = t(matrix(sapply(FIT, function(fit) fit[[2]]), nrow=ncol(MU)))

  B = t(apply(E, 2, function(y){
    b = solve(t(X) %*% X) %*% (t(X) %*% y)
  }))
  MU = X %*% t(B)

  L = lapply(FIT, function(fit) fit[[3]] - fit[[2]] %*% t(fit[[2]]) )
  SIGMA = diag(mean(diag(Reduce(`+`, L) / nrow(E))), ncol(E))

  iter = iter + 1
  print(err <- max(pB-B))
  pB = B
}


summary(m1<-lm(E~0+X))


exp(clr_coordinates(inv_ilr_coordinates(t(B))))

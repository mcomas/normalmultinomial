#2.557625 -1.529208  2.625844 -1.412314  3.460640

#2.464826 -1.409695  2.511056 -1.292941  1.803746

#2.427001 -1.373467  2.746355 -1.313165  1.818892

#2.501913 -1.432512  2.354270 -1.271511  2.150441

#2.393592 -1.399662  2.686417 -1.322034  2.111432

load('data/Pigs.rda')

X = as.matrix(Pigs)
MU = ilr_coordinates(matrix(apply(X, 2, mean), nrow=1))
SIGMA = diag(5)
library(randtoolbox)
Z = torus(1000, dim = 5, normal = TRUE)
#Z = rnormal(100, rep(0,5), diag(5))



EM_step = function(X, mu_t, sigma_t, Z){
  e = apply(X, 1, function(x){
    expect = expectedMonteCarloFixed(x, mu_t, sigma_t, Z)
    list(
      'm1' = apply(expect[[2]], 2, mean),
      'm2' = apply(expect[[3]], 1:2, mean))
  })
  MEAN  = Reduce(`+`, lapply(e, function(ee)ee[[1]])) / nrow(X)
  SIGMA = Reduce(`+`, lapply(e, function(ee)ee[[2]])) / nrow(X) - matrix(MEAN, ncol=1) %*% matrix(MEAN, nrow=1)
  list(MEAN, nearestPSD(SIGMA))
}


(PARS_T <- EM_step(X = X, mu_t = MU, sigma_t = SIGMA, Z = Z))
eigen(PARS_T[[2]])$value

for(i in 1:100){
  (PARS_T <- EM_step(X = X, mu_t = PARS_T[[1]], sigma_t = PARS_T[[2]], Z = Z))
}


E_step = function(X, mu_t, sigma_t, Z){
  t(apply(X, 1, function(x){
    expect = expectedMonteCarloFixed(x, mu_t, sigma_t, Z)
    apply(expect[[2]], 2, mean)
  }))
}
H = E_step(X = X, mu_t = PARS_T[[1]], sigma_t = PARS_T[[2]], Z = Z)
X.replaced = inv_ilr_coordinates(H)

head(X)

biplot(princomp(clr_coordinates(X.replaced)))

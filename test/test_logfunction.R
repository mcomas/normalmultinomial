library(mixpack)

set.seed(1)
N <- 10
SIZE <- 10
MU <- c(2,0)

SIGMA.ilr <- matrix(c(1,0,
                      0,1), nrow = 2)
ilrB = t(clr_coordinates(t(matrix(unlist(ilr_basis(3)), ncol=2))))
ilr_to_alr = MASS::ginv(ilrB) %*% matrix(c(1,0,-1,
                                           0,1,-1), ncol=2)
SIGMA = t(ilr_to_alr) %*% SIGMA.ilr %*% ilr_to_alr

sample = rnormalmultinomial(N, SIZE, mu = MU, sigma = SIGMA, probs = TRUE)

P = as.data.frame(sample$probs)
A = as.matrix(log(P[,1:2]/P[,3]))
X = sample$counts
X.raw = as.data.frame(X / apply(X,1,sum))

ff = sapply(1:N, function(i) mvf(a = A[i,], mu = MU, inv_sigma = solve(SIGMA), X[i,]))

ff_0 = sapply(1:N, function(i) mvf_deriv(I = 0, a = A[i,], mu = MU, inv_sigma = solve(SIGMA), X[i,]))
ff_1 = sapply(1:N, function(i) mvf_deriv(I = 1, a = A[i,], mu = MU, inv_sigma = solve(SIGMA), X[i,]))

ff_00 = sapply(1:N, function(i) mvf_deriv2(I = 0, J = 0, a = A[i,], mu = MU, inv_sigma = solve(SIGMA), X[i,]))
ff_10 = sapply(1:N, function(i) mvf_deriv2(I = 1, J = 0, a = A[i,], mu = MU, inv_sigma = solve(SIGMA), X[i,]))
ff_01 = sapply(1:N, function(i) mvf_deriv2(I = 0, J = 1, a = A[i,], mu = MU, inv_sigma = solve(SIGMA), X[i,]))
ff_11 = sapply(1:N, function(i) mvf_deriv2(I = 1, J = 1, a = A[i,], mu = MU, inv_sigma = solve(SIGMA), X[i,]))

A.max = maximize_mvf(mu = MU, inv_sigma = solve(SIGMA), X = X)
P.max = as.data.frame(cbind(exp(A.max),1) / apply(cbind(exp(A.max),1), 1, sum))

fit = normalmultinomial_fitting(X)
P.est = as.data.frame(fit[[3]])

library(ggtern)
id = LETTERS[1:N]
add_info = function(d, v) transform(d, v = v, id = id)

df = rbind(add_info(P, 'Prob.real'),
           add_info(P.max, 'Prob.max.'),
           add_info(P.est, 'Prob.est.'),
           add_info(X.raw, 'Counts'))
df$v = factor(df$v, levels = c('Prob.real', 'Counts', 'Prob.max.', 'Prob.est.'))
ggtern() +
  geom_point(data = df, aes(x=V1, y=V2, z=V3, col=id), size=3) +
  facet_wrap(~v, nrow=2) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = 'bottom')

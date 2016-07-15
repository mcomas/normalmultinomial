set.seed(1)
library(mixpack)

N = 100
MU = c(2,0)
SIZE = 10
SIGMA <- matrix(c(2,1,
                  1,2), nrow = 2)

sample = rnormalmultinomial(N, SIZE, mu = MU, sigma = SIGMA, probs = TRUE)
X = sample$counts
P = sample$probs


fit = normalmultinomial_fitting(X)
A = stepE(X, fit[[1]], fit[[2]])
P.est = cbind(exp(A), 1) / apply(cbind(exp(A),1), 1, sum)

library(ggtern)
X = as.data.frame(X)
P = as.data.frame(P)
P.est = as.data.frame(P.est)
df = bind_rows(X %>% mutate(src = 'Counts'),
               P %>% mutate(src = 'Real'),
               P.est %>% mutate(src = 'Estimated'))
ggtern() +
  geom_point(data=df, aes(x = V1, y = V2, z = V3, col = src)) +
  theme_bw() + theme(legend.position = 'none') +
  facet_wrap(~src, ncol=1)


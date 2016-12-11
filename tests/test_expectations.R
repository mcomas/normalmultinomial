library(ggtern)
library(combinat)
library(dplyr)
library(normalmultinomial)
plot_estimates = function(k){
  cat(sprintf("Sample size: %d\n", nsimplex(3,k)))
  X = t(xsimplex(3, k))
  fit_dm = dm_fit(X)
  fit_nm = nm_fit(X)
  df = bind_rows('count' = data.frame(X),
                 'dm' = data.frame(fit_dm$expected),
                 'nm' = setNames(data.frame(fit_nm$expected), c('X1','X2','X3')), .id='meth')
  ggtern(data = df) +
    geom_mask() +
    geom_point(aes(x=X1, y=X2, z=X3), size=1) +
    facet_wrap(~meth, nrow=1) +
    theme_classic() +
    ggtitle(sprintf('(3,%d)-simplex-lattice', k))
}

plot_estimates(1)
plot_estimates(2)
plot_estimates(3)
plot_estimates(4)
plot_estimates(5)
plot_estimates(10)
plot_estimates(15)
plot_estimates(20)
plot_estimates(30)


X = matrix(c(3,2,1,
             3,1,2,
             2,1,3,
             2,3,1,
             1,2,3,
             1,3,2), ncol = 3, byrow = T)
fit_dm = dm_fit(X)
fit_nm = nm_fit(X, verbose = TRUE, eps = 10^-4)


X = matrix(c(0,2,1,
             0,1,2,
             2,1,0,
             2,0,1,
             1,2,0,
             1,0,2), ncol = 3, byrow = T)
fit_dm = dm_fit(X)
fit_nm = nm_fit(X, verbose = TRUE, eps = 10^-4)

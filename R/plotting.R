library(Rcpp)
sourceCpp('src/fitting.cpp')

if(!exists('PAR')){
  PAR = list(seed = 150, min_int = 2.5, max_int = 3.5, Mu = 4, Sigma = matrix(3), Size = 1000)
  PAR = list(seed = 150, min_int = -100, max_int = 25, Mu = -4, Sigma = matrix(300), Size = 1000)
}

set.seed(PAR$seed)
Mu = PAR$Mu
Sigma = PAR$Sigma
Size = PAR$Size

(X <- as.numeric( rnormalmultinomial(mu = Mu, size = Size, sigma = Sigma) ) )
#X = as.numeric( rnormalmultinomial(mu = Mu, size = Size, sigma = Sigma) )
#X = as.numeric( rnormalmultinomial(mu = Mu, size = Size, sigma = Sigma) )
#X = as.numeric( rnormalmultinomial(mu = Mu, size = Size, sigma = Sigma) )
#X = as.numeric( rnormalmultinomial(mu = Mu, size = Size, sigma = Sigma) )

step = 0.01
A = seq(PAR$min_int, PAR$max_int, step)
cLog = sapply(A, mvf, mu = Mu, inv_sigma = 1/Sigma, x = X)

df = data.frame(A = as.numeric(A), LogLik = as.numeric(cLog))
df$dens = exp(df$LogLik) / sum(exp(df$LogLik) * step)
df$deriv = sapply(df$A, function(a) mvf_deriv(0, a, mu = Mu, inv_sigma = 1/Sigma, x = X))
df$hess = sapply(df$A, function(a) mvf_deriv2(0, 0, a, mu = Mu, inv_sigma = 1/Sigma, x = X))

(E_A <- sum(df$A*df$dens * step) )
(E_AA <- sum(df$A*df$A*df$dens * step) )
(E_kappa <- sum(log(1+exp(df$A))*df$dens * step) )


(Max_A <- as.numeric(Mstep(A=matrix(E_A), mu = Mu, inv_sigma = solve(Sigma), X = matrix(X, nrow=1))))


library(ggplot2)
library(dplyr)
library(reshape2)

f = function(a) approx(x=df$A, y=df$dens, xout = a)$y

d = df %>% mutate(logdens = log(dens), expLogLik = exp(LogLik)) %>% select(A, dens, logdens, expLogLik, LogLik,  deriv, hess) %>% melt(id.vars = c('A'))
ggplot(data=d) + 
  geom_segment(y=0, yend=f(E_A), x=E_A, xend=E_A, col='red') + 
  geom_segment(y=0, yend=f(Max_A), x=Max_A, xend=Max_A, col='blue') +
  geom_line(aes(x=A, y=value)) +
  facet_grid(variable~., scales = 'free')



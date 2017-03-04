library(normalmultinomial)
library(randtoolbox)

x = c(1,0)
MU = 0
SIGMA = matrix(1)

NSIM = 1000
MU_EXP = 0.5135884080
SIGMA_EXP = 1

set.seed(1)
fit0 <- expectedMonteCarlo(x, MU, SIGMA,
                           Z = matrix(rnorm(NSIM/2), ncol=1),
                           mu_exp = MU_EXP)

set.seed(1)
fit1 <- expectedMonteCarlo_withoutAV(x, MU, SIGMA,
                           Z = matrix(rnorm(NSIM), ncol=1),
                           mu_exp = MU_EXP)

set.seed(1)
fit2 <- expectedMetropolis(x, MU, SIGMA,
                           Z = matrix(rnorm(NSIM), ncol=1), mu_exp = MU_EXP)

m1.pseudo = fit0[[2]]
m1.withoutAV = fit1[[2]]
m1.metropolis = fit2[[2]]

m2.pseudo = c(fit0[[3]])
m2.withoutAV = c(fit1[[3]])
m2.metropolis = c(fit2[[3]])

pdf('figures/montecarlo_comparison.pdf', width = 8.3, height = 3.8)
par(mfrow=c(1,2))
plot(rep(cumsum(m1.pseudo)/seq_along(m1.pseudo), each=2),
     type='l', xlab='Iterations', ylab = 'Estimation',
     main='First moment', ylim =  MU_EXP + c(-1,1) * 0.1)
lines(cumsum(m1.withoutAV)/seq_along(m1.withoutAV), col='blue')
lines(cumsum(m1.metropolis)/seq_along(m1.metropolis), col='green')
abline(h=MU_EXP, col='red')

legend('topright', c('Pseudo generation', 'Without AV', 'Metropolis'),
       col=c('black', 'blue', 'green'), bty = 'n', lty=1, cex = 0.8)


plot(rep(cumsum(m2.pseudo)/seq_along(m2.pseudo), each=2),
     type='l', xlab='Iterations', ylab = 'Estimation',
     main='Second moment', ylim = SIGMA_EXP + c(-1,1) * 0.15)
lines(cumsum(m2.withoutAV)/seq_along(m2.withoutAV), col='blue')
lines(cumsum(m2.metropolis)/seq_along(m2.metropolis), col='green')

abline(h=1, col='red')
legend('topright', c('Pseudo generation', 'Without AV', 'Metropolis'),
       col=c('black', 'blue', 'green'), bty = 'n', lty=1, cex = 0.8)
dev.off()


res = list()
set.seed(1)
res[['Pseudo generation']] = replicate(500, sapply(
  expectedMonteCarlo(x, MU, SIGMA,
                     Z = matrix(rnorm(NSIM)/2, ncol=1),
                     mu_exp = MU_EXP),
  mean))
set.seed(1)
res[['Without AV']] = replicate(500, sapply(
  expectedMonteCarlo_withoutAV(x, MU, SIGMA,
                               Z = matrix(rnorm(NSIM), ncol=1),
                               mu_exp = MU_EXP),
  mean))
set.seed(1)
res[['Metropolis']] = replicate(500, sapply(
  expectedMetropolis(x, MU, SIGMA,
                     Z = matrix(rnorm(NSIM), ncol=1), mu_exp = MU_EXP),
  mean))

tab2.mean = t(sapply(res, function(r) setNames(apply(r,1,mean)[2:3], c('M1','M2'))))
tab2.sd = t(sapply(res, function(r) setNames(apply(r,1,sd)[2:3], c('M1','M2'))))

tab.out = data.frame(
  'Method' = gsub('_','\\_',rownames(tab2.mean), fixed=T),
  'First moment' = sprintf("%8.7f (%6.5f)", tab2.mean[,1], tab2.sd[,1]),
  'Second momenth' = sprintf("%8.7f (%6.5f)", tab2.mean[,2], tab2.sd[,2]))
tab.out

library(Hmisc)
sink(file='tex/convergence01.tex')
cat(latexTabular(tab.out,
                 headings = sprintf("\\textbf{%s}", names(tab.out)),
                 hline = 1,
                 align = 'r | r r',
                 helvetica = F,
                 translate = F))
sink()

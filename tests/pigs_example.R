load('data/Pigs.Posture.RData')

library(zCompositions)
library(normalmultinomial)
source('tests/biplots.R')
## Load data

load("data/Pigs.Posture.RData")

comp.pos <- dat.pos[,3:8]
zPatterns(comp.pos,0) # 16.09% overall zeros, 6 patterns, some complexity
dev.off()
X = as.matrix(comp.pos)

## Repl using normal-multinomial
set.seed(1)
library(parallel)
cl = makeCluster(detectCores()-2)
fit = nm_fit(X, verbose = T, parallel.cluster = cl)
P.nm <- fit$expected
colnames(P.nm) = colnames(X)
stopCluster(cl)

## Repl using dirichlet-multinomial
fit.dm = dm_fit(X)
P.dm <- fit.dm$expected
colnames(P.dm) = colnames(X)

## Repl using GBM
P.gbm <- cmultRepl(comp.pos)

## Repl using GBM
P.user <- cmultRepl(as.matrix(comp.pos),
                    method = 'user',
                    t=matrix(fit.dm$alpha/sum(fit.dm$alpha),
                             ncol=ncol(comp.pos),
                             nrow=nrow(comp.pos),
                             byrow = T),
                    s=rep(sum(fit.dm$alpha),nrow(comp.pos)))


with0 = ifelse(X[,5] == 0, 'with 0', 'without 0')
g = paste(dat.pos$Treatment, with0, sep='-')

## Compare biplots
library(ggbiplot)
# ordinary biplot, no replacement
ggbiplot(princomp(comp.pos),groups=g) +
  coord_fixed(xlim = c(-2.5,2.5), ylim=c(-2.5,2.5)) +
  ggtitle('Ordinary biplot')

# ordinary corr biplot, no replacement
ggbiplot(princomp(scale(comp.pos,center=FALSE)),groups=g) +
  coord_fixed(xlim = c(-2.5,2.5), ylim=c(-2.5,2.5)) +
  ggtitle('Ordinary corr biplot \n(no replacement)')

ggbiplot(princomp(clr_coordinates(P.gbm)),groups=g) +
  coord_fixed(xlim = c(-2.5,2.5), ylim=c(-2.5,2.5)) +
  ggtitle('Geometric Bayesian\n(multiplicative replacement)')

ggbiplot(princomp(clr_coordinates(P.dm)),groups=g) +
  coord_fixed(xlim = c(-2.5,2.5), ylim=c(-2.5,2.5)) +
  ggtitle('Dirichlet-Multinomial')

ggbiplot(princomp(clr_coordinates(P.user)),groups=g) +
  coord_fixed(xlim = c(-2.5,2.5), ylim=c(-2.5,2.5)) +
  ggtitle('Dirichlet-Multinomial \n(multiplicative replacement)')

ggbiplot(princomp(clr_coordinates(P.nm)),groups=g) +
  coord_fixed(xlim = c(-2.5,2.5), ylim=c(-2.5,2.5)) +
  ggtitle('Normal-Multinomial')

# Parece que los CoDa ayudan a distinguir mejor los grupos.


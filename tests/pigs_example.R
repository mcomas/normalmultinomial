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
fit = nm_fit(X)
P1 <- nm_expected(X, fit$mu, fit$sigma)
colnames(P1) = colnames(X)

## Repl using GBM
P2 <- cmultRepl(comp.pos)

with0 = ifelse(X[,5] == 0, 'with 0', 'without 0')
g = paste(dat.pos$Treatment, with0, sep='-')

## Compare biplots
library(ggbiplot)
ggbiplot(princomp(comp.pos),groups=g) # ordinary biplot, no replacement
ggbiplot(princomp(scale(comp.pos,center=FALSE)),groups=g) # ordinary corr biplot, no replacement
ggbiplot(princomp(clr_coordinates(P1)),groups=g)
ggbiplot(princomp(clr_coordinates(P2)),groups=g)

# Parece que los CoDa ayudan a distinguir mejor los grupos.


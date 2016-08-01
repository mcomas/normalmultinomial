library(normalmultinomial)
load('data/Pigs.rda')
Pigs = as.matrix(Pigs)

E = nm_fit(Pigs)
max(E[[2]])

E = nm_fit(Pigs[,c(1,3,5)])
max(E[[2]])

nm_fit(Pigs)[[1]][,2:7]


X
nm_fit(X)

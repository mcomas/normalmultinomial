library("MASS")
data(menarche)

Y = with(menarche, cbind(Menarche, Total-Menarche))

mod = glm(Y~1, family = binomial(logit))


fit = nm_fit(Y, eps = 0.00000001, verbose = TRUE)

plot(1:nrow(Y), Y[,1]/(Y[,1]+Y[,2]))
points(1:nrow(Y), mod$fitted.values)
points(1:nrow(Y), fit$expected[,1], col='red')
points(1:nrow(Y), rep(inv_ilr_coordinates(fit$mu)[1], nrow(Y)), col='red')

y = c(ilr_coordinates(fit$eYpected))
summary(m1<-glm(Y~Age, family = binomial(logit), data=menarche))
summary(m2<-lm(y~Age, data=menarche))


###
# Parliament
data("parliament2015")
X = parliament2015[, c('cup', 'cs', "birth.cat", "pop")]
X$birth.cat = X$birth.cat / X$pop
X = X[order(X$birth.cat),]

mod = glm(cbind(cup,cs)~birth.cat, family = binomial(logit), data=X)
fit = nm_fit(as.matrix(X[,c('cup','cs')]), eps = 0.00001, verbose = TRUE)

plot(1:nrow(X), X[,1]/(X[,1]+X[,2]))
points(1:nrow(X), mod$fitted.values)
points(1:nrow(X), fit$expected[,1], col='red')
points(1:nrow(X), rep(inv_ilr_coordinates(fit$mu)[1], nrow(X)), col='red')


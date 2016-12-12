library("MASS")
data(menarche)

X = with(menarche, cbind(Menarche, Total-Menarche))

mod = glm(X~1, family = binomial(logit))


fit = nm_fit(X, eps = 0.00000001, verbose = TRUE)

plot(1:nrow(X), X[,1]/(X[,1]+X[,2]))
points(1:nrow(X), mod$fitted.values)
points(1:nrow(X), fit$expected[,1], col='red')
points(1:nrow(X), rep(inv_ilr_coordinates(fit$mu)[1], nrow(X)), col='red')

y = c(ilr_coordinates(fit$expected))
summary(glm(X~Age, family = binomial(logit), data=menarche))
summary(lm(y~Age, data=menarche))

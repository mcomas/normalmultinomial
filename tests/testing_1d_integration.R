
f = function(h, x, mu, sigma, m, center){
  logf = dnorm(h, mu, sigma, log=TRUE) + dbinom(x[1], sum(x), exp(h * sqrt(0.5)) /  (exp(h * sqrt(0.5)) + exp(-h * sqrt(0.5))), log = TRUE)
  (h-center)^m * exp(logf)
}


Px = function(X, mu, sigma, m, center = rep(0, nrow(X))){
  sapply(1:nrow(X), function(i){
    x = X[i,]
    xs = seq(-100, 100, length=1000)
    rank = f(xs, x = x, mu = mu, sigma=sigma, m=m, center=center[i])
    rank.v = range(xs[rank!=0])
    linf = max(xs[xs < rank.v[1]])
    lsup = min(xs[xs > rank.v[2]])
    # # Version 1
    # h = seq(linf, lsup, length.out = 100000)
    # sum(f(h, x = x, mu = mu, sigma=sigma, m=1)*(diff(range(h))/length(h)))
    # Version 2
    integrate(f, linf, lsup, x = x, mu = mu, sigma=sigma, m = m, center=center[i])$value
  })
}



(mu <- ilr_coordinates( colMeans(X/rowSums(X)) ))
sigma = 1

for(i in 1:10){
  m0 = Px(X, mu = mu, sigma=sigma, m = 0)
  m1 = Px(X, mu = mu, sigma=sigma, m = 1) / m0
  m2 = Px(X, mu = mu, sigma=sigma, m = 2, center = m1) / m0

  print(mu <- mean(m1))
  print(sigma <- sqrt(mean(m2)))
}

for(i in 1:nrow(X)){
  print(i)
  x = X[i,]
  xs = seq(-100, 100, length=10000)
  rank = f(xs, x = x, mu = mu, sigma=sigma, m=m)
  rank.v = range(xs[rank!=0])
  linf = max(xs[xs < rank.v[1]])
  lsup = min(xs[xs > rank.v[2]])
  # # Version 1
  # h = seq(linf, lsup, length.out = 100000)
  # sum(f(h, x = x, mu = mu, sigma=sigma, m=1)*(diff(range(h))/length(h)))
  # Version 2
  integrate(f, linf, lsup, x = x, mu = mu, sigma=sigma, m=m)$value
  i = i +1
}

integrate(f, -Inf, Inf, x = x, mu = mu, sigma=sigma, m=1)

h = seq(0, 300, length.out = 10000)
sum(f(h, x = x, mu = mu, sigma=sigma, m=1)*(diff(range(h))/length(h)))

xs = seq(-100, 100, length=100)
rank = f(xs, x = x, mu = mu, sigma=sigma, m=1)
rank.v = range(xs[rank!=0])
linf = max(xs[xs < rank.v[1]])
lsup = min(xs[xs > rank.v[2]])
# Version 1
h = seq(linf, lsup, length.out = 100000)
sum(f(h, x = x, mu = mu, sigma=sigma, m=1)*(diff(range(h))/length(h)))
# Version 2
integrate(f, linf, lsup, x = x, mu = mu, sigma=sigma, m=1)$value


integrate(f, lower = ilr_coordinates(x)-sdev*sigma, upper =ilr_coordinates(x)+sdev*sigma, x=x, m=m)


  adaptIntegrate(f, lowerLimit = mu-sdev*sigma, upperLimit = mu+sdev*sigma, x=x, m=m, maxEval = 1000)

library(microbenchmark)

X = diag(3)
MU = matrix(c(0,0), 1)
SIGMA = diag(2)

df_x_1(X, MU, SIGMA, nsim = 100000)
df_x_2(X, MU, SIGMA, nsim = 100000)
df_x_3(X, MU, SIGMA, nsim = 100000)

microbenchmark(
  all = df_x_1(X, MU, SIGMA, nsim = 1000),
  ind = df_x_2(X, MU, SIGMA, nsim = 1000),
  max = df_x_3(X, MU, SIGMA, nsim = 1000) )

N = 1000
(E <- expected1(X, MU, SIGMA, nsim = N))[[1]]
round(inv_ilr_coordinates(E[[1]][,2:3]),3)
sqrt(E[[2]]/N)

(E <- expected_initial(X, MU, SIGMA, se_eps = 0.0001))[[1]]
round(inv_ilr_coordinates(E[[1]][,2:3]),3)


(E <- expected1(X, MU, SIGMA, nsim = 1000))
round(inv_ilr_coordinates(E[[1]][,2:3]),3)
(E <- expected2(X, MU, SIGMA, se_eps = 0.001))
round(inv_ilr_coordinates(E[[1]][,2:3]),3)
(E <- expected3(X, MU, SIGMA, se_eps = 0.001))
round(inv_ilr_coordinates(E[[1]][,2:3]),3)


M0 = t(E[[1]][,2:3])
S0 = array(t(E[[1]][,c(4,5,5,6)]), dim = c(2,2,3)) -
  array(apply(M0, 2, function(col) col %*% t(col)), dim=c(2,2,3))
(E <- expected4(X, MU, SIGMA, M0, S0, se_eps = 0.001))
round(inv_ilr_coordinates(E[[1]][,2:3]),3)

(E <- expected5(X, MU, SIGMA, M0, S0, se_eps = 0.001))
round(inv_ilr_coordinates(E[[1]][,2:3]),3)

(E <- expected6(X, MU, SIGMA, M0, S0, se_eps = 0.001))
round(inv_ilr_coordinates(E[[1]][,2:3]),3)




#####
####
###
### Extrange thinks happen
X = matrix(c(3,2,1,
             3,1,2,
             2,1,3,
             2,3,1,
             1,2,3,
             1,3,2), ncol = 3, byrow = T)
(E <- expected1(X, MU, SIGMA, nsim = 1000))
round(inv_ilr_coordinates(E[[1]][,2:3]),3)
b = matrix(c(1,2,3), ncol = 1)
ALR = t(t(log(b[1:2,]/b[3,])))
ILR = t(ilr_coordinates(t(b)))

solve(ilr_to_alr(3)) %*% ALR

########
(E <- expected3(matrix(c(1,0,0,0,1,0,0,0,1), ncol=3), MU, SIGMA, se_eps = 0.001))

M0 = t(E[[1]][,2:3])
S0 = array(t(E[[1]][,c(4,5,5,6)]), dim = c(2,2,3)) -
  array(apply(M0, 2, function(col) col %*% t(col)), dim=c(2,2,3))
(E <- expected4(matrix(c(1,0,0,0,1,0,0,0,1), ncol=3, byrow = T), MU, SIGMA,
                M0, S0, se_eps = 0.001))
round(inv_ilr_coordinates(E[[1]][,2:3]),3)

(E <- expected5(matrix(c(1,0,0,0,1,0,0,0,1), ncol=3, byrow = T), MU, SIGMA,
                M0, S0, se_eps = 0.001))
round(inv_ilr_coordinates(E[[1]][,2:3]),3)

(E <- expected6(matrix(c(1,0,0,0,1,0,0,0,1), ncol=3, byrow = T), MU, SIGMA,
                M0, S0, se_eps = 0.0001))
round(inv_ilr_coordinates(E[[1]][,2:3]),3)

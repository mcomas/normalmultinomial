library(mixpack)

set.seed(1)
MU <- c(0,0)
SIGMA.ilr <- matrix(c(1,0,
                      0,1), nrow = 2)
ilrB = t(clr_coordinates(t(matrix(unlist(ilr_basis(3)), ncol=2))))
ilr_to_alr = MASS::ginv(ilrB) %*% matrix(c(1,0,-1,
                                           0,1,-1), ncol=2)
SIGMA = t(ilr_to_alr) %*% SIGMA.ilr %*% ilr_to_alr


choose_starting(x = c(0,20,0), mu = MU, inv_sigma = solve(SIGMA), prop = 0.06)

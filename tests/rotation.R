
D = 3
x = matrix(0, nrow = D, ncol = D+1)
for(i in 1:D){
  s = sum(x[1:i,i]^2)
  x[i,i] = sqrt ( 1 - s )
  for(j in 1+i:D){
    s = 0
    for(ii in 1:i){
      s = s + x[ii,i] * x[ii,j]
    }
    x[i,j] = ( -1/D - s ) / x[i,i]
  }
}

x

t(ilr_basis(4)/ilr_basis(4)[4,3])[D:1,(D+1):1]

v1 = c(1,2,3,4)
v2 = c(3,5,2,4)



# 1. Let v3 = v1 X v2.
# 2. Let Q2 be an orthogonal matrix that performs a rotation about v1
# by an arbitrary angle so that Q2 * v1 = v1. This means that v1 is an
# eigenvector of Q2.
# 3. Let Q1 be an orthogonal matrix that performs a rotation about v3
# by the smallest angle between v1 and v2 so that Q1 * v1 = v2.
#
# We can see that Q1 * (Q2 * v1) = v2.  Now define M = Q1 * Q2 so that
#
# M * v1 = v2
#
# Since there are an infinite number of options for Q2, there are an
# infinite number of options for M.  The simplest option is when Q2 = I.
#
# Now you just need to know how to construct an orthogonal matrix for a
# rotation about an axis by a specified angle.  For this, you need
# Rodrigues' Rotation Formula.  Here are two articles that explain it:

#LU function
library(Matrix)
lu_function <- function (X, y) {
  XX <- crossprod(X)
  XX.lu <- expand(lu(XX))
  L <- XX.lu$L
  U <- XX.lu$U
  P <- XX.lu$P
  bs <- forwardsolve(L, crossprod(P, crossprod(X, y)))
  betahat <- backsolve(U, bs)
  
  #sigma^2 hat
  #sigma2hat=crossprod(y - (X) %*% betahat) / (dim(X)[1] - dim(X)[2]-1)
  
  #return answer
  list(coefficients = betahat)
}


#solve (In-P')x = 0n by using SVD decomposition
svd_function <- function (X, y) {
  
  #  Get svd for X matrix
  svd_mat=svd(X)
  svd_sing=svd_mat$d
  U_mat <- svd_mat$u
  V_mat <- svd_mat$v
  
  # Solve x = V*(1/d)*u'*y to get the all beta values
  beta <- as.vector(V_mat%*%diag(1/(svd_sing))%*%t(U_mat)%*%y)
  
  # return answer
  list(coefficients = beta)
}
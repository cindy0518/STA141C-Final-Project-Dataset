#QR 

qr_function <- function (X, y) {
  #  Get Q and R matrix for X matrix
  qr_mat=qr(X)
  Q_mat <- qr.Q(qr_mat)
  R_mat <- qr.R(qr_mat)
  
  # Solve Rb = Q'y to get the all beta values
  betahat <- as.vector(backsolve(R_mat, crossprod(Q_mat, y)))
  
  # sigma^2 hat
  #sigma2hat=crossprod(y - (X) %*% betahat) / (dim(X)[1] - dim(X)[2]-1)
  
  # return answer
  list(coefficients = betahat)
}

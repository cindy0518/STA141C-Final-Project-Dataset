Jac_iter=function(A, maxiter = 1000,tolx = 1e-3) {
  n <- dim(A)[1]
  A=as.matrix(A)
  A_sparse <- Matrix(A, sparse = TRUE)
  outd_s <- rowSums(A_sparse)
  outdinv_s <- ifelse(outd_s > 0, 0.85 / outd_s, 0)
  z <- ifelse(outd_s > 0, 0.15, 1) / n
  # G is the R^+ A matrix
  G <- outdinv_s * A_sparse
  # set starting point
  x = rep(1 / n, n)
  # obtain diagonal of the linear system
  aii <- 1 - diag(G) - z
  # Jacobi iterations
  for (iter in 1:maxiter) {
    xold <- x
    # this is (I - P^T)x
    x <- x - x %*% G - sum(x * z)
    # this is Jacobi update
    x <- xold - x / aii
    if (sum(abs(x - xold)) < tolx)
      break
  }
  if (iter == maxiter) {
    warning(paste("fail to converge in", maxiter, " iterations!"))
  } else {
    print(paste("converged in", iter, "iterations"))
  }
  pgrank <- abs(x / sum(x))
  return(pgrank)
}
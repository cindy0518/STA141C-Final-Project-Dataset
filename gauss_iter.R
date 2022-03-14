gauss_literation =function(A, telep = 0.85,  maxiter = 1000,
         tolx = 1e-3) {
  n<- dim(A)[1]
  ##REFERENCE STA 141C HW-2 Notes:
  A=as.matrix(A)
  A_sparse <- Matrix(A, sparse = TRUE)
  outd_s <- rowSums(A_sparse)
  outdinv_s <- ifelse(outd_s > 0, telep / outd_s, 0)
  z <- ifelse(outd_s > 0, 1 - telep, 1) / n
  # G is the R^+ A matrix
  G <- outdinv_s * A_sparse
  #set the initial vector 
  x0 = rep(1 / n, n)
  #b from Ax = b:
  bi = x0
  #Diagonal of the Matrix:
  aii <- 1 - diag(G) - z
  k = 1
  while(k < maxiter){
    x = x0
    Gauss_dec <- bi - x0 %*% G - sum(x * z)
    x <- bi - Gauss_dec / aii
    if (any(abs(x - x0) < tolx)){
      gs_pgrank<- x
      break
    }
    k = k +1
    
    for (i in 1: maxiter){
      x0 = x
    }
  }
  #normalize
  gs_pgrank <- abs(x0 / sum(x0))
  #output
  return(gs_pgrank)
}

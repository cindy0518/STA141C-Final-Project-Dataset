literation = function(A) {
  
  matrixdim = dim(A)[1]
  
  #create the list to input the q vector
  q = list()
  #initial guess of q^0 
  q[[1]] = matrix(c(rep(1 / matrixdim, matrixdim)), nrow = matrixdim)
  
  #iterative q^k and q^(k-1) by using q^k=Aq^(k-1)/norm(Aq^(k-1))
  for (i in 1:10000) {
    #here we using i+1 instead of i because the first position in q is initial value
    q[[i + 1]] = A %*% q[[i]] / norm(A %*% q[[i]], type = c("2"))
    #length for q^k and q^(k-1)
    length_x = norm(q[[i + 1]] - q[[i]], type = c("2"))
    
    #condition for return value of q^k, lambda, iteration times 
    #when length_x<= tolerance 10^-6
    if (10 ^ -3 >= length_x) {
      eigenvector_q=q[[i+1]]
      lamdba = t(q[[i + 1]]) %*% A %*% q[[i]]
      return(list(times=length(q) - 1,lamdba=lamdba,eigenvector_q=eigenvector_q/sum(eigenvector_q)))
    }
  }
}
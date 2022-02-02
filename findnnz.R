findnnz <- function(x,r,tiny=1e-16){
  n <- length(x)
  normx1 <- norm(x)
  #print(normx1)
  
  if (min(normx1, norm(x, type='I')) <= tiny){
    k <- 0
    xnew <- Matrix(0, n, 1)
    return(list(k=k, xnew=xnew))
  }
  
  absx <- abs(x[order(abs(x), decreasing=TRUE)])
  idx <- order(abs(x), decreasing=TRUE)
  tmpidx <- which(cumsum(absx) > r*normx1)
  
  #print(absx)
  #print(idx)
  
  k <- tmpidx[1]
  #print(k)
  xnew <- Matrix(0, n, 1)
  idxnew <- idx[1:k]
  xnew[idxnew] <- x[idxnew]
  return(list(k=k, xnew=xnew))
}
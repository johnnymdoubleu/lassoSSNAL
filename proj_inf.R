proj_inf <- function(x, lambda){
  if (all(lambda < 0 )){
    print("error")
  }
  else if(sum(lambda) == 0){
    y <- 0 * x
  }
  else{
    y <- pmax(-lambda, pmin(x,lambda))
  }
  rr <- (y==x)
  rr <- Matrix(rr) * 1
  return(list(y=y, rr=rr))
}
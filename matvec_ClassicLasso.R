matvec_Classic_Lasso <- function(y, par, AP){
  #temp <- t(AP) %*% y
  #Ay <- y + par$sigma * (AP %*% temp)
  #print('--')
  #print(dim(AP))
  
  temp <- eigenMatMult(t(AP),y,4)
  Ay <- y + par$sigma * eigenMatMult(AP,temp,4)
  return(Ay)
}

matvec_Classic_Lasso_Amap <- function(y, par, Ainput) {
  temp <- Ainput$ATmpa[y]
  temp <- (1-par$rr) * temp
  Ay <- y + par$sigma * Ainput$Amap[temp]
  return(Ay)
}
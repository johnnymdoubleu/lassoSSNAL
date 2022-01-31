matvec_Classic_Lasso <- function(y, par, AP){
  temp <- t(AP) %*% y
  Ay <- y + par$sigma * (AP %*% temp)
}

matvec_Classic_Lasso_Amap <- function(y, par, Ainput) {
  temp <- Ainput$ATmpa[y]
  temp <- (1-par$rr) * temp
  Ay <- y + par$sigma * Ainput$Amap[temp]
}
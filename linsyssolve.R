# library to call mldivide()
library(pracma)

# additional function declared to 
speye <- function(N){
  return(as.matrix(sparseMatrix(i=(1:N),j=(1:N), x=rep(1,N) ) ) )
}

linsyssolve <- function(Ainput, rhs, par){
  m <- length(rhs)
  pp <- as.logical(!par$rr)
  #Ayes <- exists(A)
  Ayes<-1
  solve <- "d_pcg"
  dn <- 10000
  sp <- sum(pp)
  
  #return(sp)
  
  if (m<=dn & Ayes){
    if (m <= 1000){
      solver <- 'd_direct'
    }
    else if(sp <= max(0.01 * par$n, dn)){
      solver <- "d_direct"
    }
  }
  if (sp <= 0.7*m & Ayes & sp <= dn){
    solver <- "p_direct"
  }
  if ((m>5e3 & sp >= 200)| (m>2000 & sp > 800) | (m > 100 & sp > 1e4)){
    solver <- "d_pcg"
  }
  
  print(solver)
  
  if (solver == "d_pcg") {
    if (Ayes){
      #AP <- Ainput$A[,pp]
      AP <- Ainput[,pp]
      #print("solverapsum")
      #print(sum(AP))
      if (FALSE){
        tmp <- sum(AP * t(AP))
        par$precond <- 1
        par$invdiagM <- 1 / (1 + par$sigma*tmp)
      }
      #c(xi,psq, resnrm, solve_ok) <- psqmry('matvec_classicLasso_Amap',AP, rhs, par)
      #return(psqmry('matvec_Classic_Lasso',AP, rhs, par))
      psqmryout <- psqmry('matvec_Classic_Lasso', AP, rhs, par)
      output <- list(xi = psqmryout$x,
                     resnrm = psqmryout$resnrm,
                     solve_ok = psqmryout$solve_ok)
      return(output)
    }
    else {
      c(xi,psq, resnrm, solve_ok) <- psqmry('matvec_classicLasso_Amap',Ainput, rhs, par)
    }
  }
  else if(solver == "d_direct"){
    # AP <- Ainput$A[,pp]
    AP <- Ainput[,pp]
    sigAPAt <- par$sigma*(AP%*%t(AP))
    # sigAPAt <- par$sigma * eigenTransMapMatMult(AP, 4)
    if (m <= 1500){
      M <- diag(m) + sigAPAt
      xi <- mldivide(M, rhs) # same as backslash operator in Matlab
    }
    else {
      M <- speye(m) + sigAPAt
      L <- chol(M)
      xi <- solve(L,rhs)
    }
    resnrm <- 0
    solve_ok <- 1
  }
  else if(solver == "p_direct"){
    #AP <- Ainput$A[,pp]
    AP <- Ainput[,pp]
    #print(dim(A))
    #print(pp)
    #print(par$rr)
    #AP <- A[,(pp+1)]
    APT <- t(AP)
    # rhstmp <- APT %*% rhs
    rhstmp <- eigenMapMatMult(APT, rhs, 4)

    
    # PAtAP <- APT %*% AP
    PAtAP <- eigenMapMatMult(APT, AP,4)
    if (sp <= 1500) {
      M <- diag(sp)/par$sigma + PAtAP
      if(0 %in% dim(M)) {
        tmp <- matrix(NA,0,1)
        #return(dim(tmp))
      } 
      else {
        tmp <- mldivide(M, rhstmp)
      }
    }
    else {
      M <- speye(sp)/par$sigma + PAtAP
      L <- chol(M)
      tmp <- solve(L, rhstmp)
    }
    resnrm <-0
    solve_ok <- 1
    # xi <- rhs - eigenVecMatMult(AP, tmp, 4)
    xi <- rhs - AP %*% tmp
    
    #print("here")
  }
  
  output <- list(xi = xi,
                 resnrm = resnrm,
                 solve_ok = solve_ok)
  return(output)
}
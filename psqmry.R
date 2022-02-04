psqmry <- function(matvec, A, b, par, x0=NULL, Ax0=NULL){
  
  N <- length(b)
  maxit <- max(5000, sqrt(N))
  tol <- 1e-6 * norm(b, "2")
  stagnate_check <- 20
  miniter <- 0
  
  #print(tol)
  
  # nargin <- length(as.list(match.call()))-1
  # if(nargin < 5){
  #   x0 <- matrix(0, N)
  # }
  # I am using !exists() here as nargin is not available in R
  
  if (is.null(x0)){x0 <- matrix(0, N)}
  if ("maxit" %in% names(par)){maxit <- par$maxit}
  if ("tol" %in% names(par)){tol <- par$tol}
  if ("stagnate_check_psqmr" %in% names(par)){stagnate_check <- par$stagnate_check_psqmr}
  if ("minitpsqmr" %in% names(par)){miniter <- par$minitpsqmr}
  
  solve_ok <- 1
  
  printlevel <- 0

  x <- x0
  
  #print(norm(x,"2"))
  
  if (norm(x,"2") > 0){
    if (!is.null(Ax0)){
      Ax0 <- feval(matvec, x0, par$sigma, A)
    }
    Aq <- Ax0
  }
  else {
    Aq <- matrix(0, N, 1)
  }
  r <- b - Aq
  err <- norm(r,"2")
  resnrm <- err
  minres <- err
  
  #print(err)
  ##
  q <- precondfun(par,r)
  
  #print(sum(q))
  
  tau_old <- norm(q,"2")
  rho_old <- as.numeric(t(r) %*% q) #might use eigenTransMatMult()
  
  #print(tau_old)
  #print(rho_old)
  
  theta_old <- 0
  d <- matrix(0,N,1)
  res <- r
  Ad <- matrix(0,N,1)
  
  ##
  ## Main Loop
  ##
  tiny <- 1e-30
  for (iter in 1:maxit){
    #print(dim(A))
    
    Aq <- feval(matvec, q, par$sigma, A)
    sigma <- as.numeric(t(q) %*% Aq) #might use eigenTransMatMult()
    #print(matvec)
    #print(sigma)
    #print(sum(Aq))
    
    if (abs(sigma) < tiny){
      solve_ok <- 2
      if (printlevel) cat('s1')
      break
    }else{
      alpha <- rho_old/sigma
      r <- r - alpha*Aq
    }
    u <- precondfun(par,r)
    
    ##
    theta <- norm(u, "2")/tau_old
    c <- 1/sqrt(1+theta^2)
    tau <- tau_old*theta*c
    gam <- (c^2*theta_old^2)
    eta <- c^2*alpha
    d <- gam*d + eta*q
    x <- x + d
    ##-------------stopping conditions---------
    Ad <- gam*Ad + eta*Aq
    res <- res - Ad
    err <- norm(res,"2")
    resnrm <- c(resnrm,err)
    if (err < minres) minres <- err
    if ((err < tol) && (iter > miniter) && (t(b) %*% x > 0)) break
    if ((iter > stagnate_check) && (iter > 10)){
      ratio = resnrm[(iter-9):(iter+1)]/resnrm[(iter-10):iter]
      if ((min(ratio) > 0.997) && (max(ratio) < 1.003)){
        if (printlevel) cat('s')
        solve_ok <- -1
        break
      }
    }
    ##--------------------------------------------
    if (abs(rho_old) < tiny){
      solve_ok <- 2
      cat('s2')
      break
    }else{
      rho <- as.numeric(t(r) %*% u) #possible use of eigenTransMatMult()
      beta <- rho/rho_old
      q <- u + beta*q
    }
    rho_old <- rho
    tau_old <- tau
    theta_old <- theta
  }
  if (iter == maxit) solve_ok <- -2
  if (solve_ok != -1){
    if (printlevel) cat(' ')
  }
  
  return (list("x"=x,"resnrm"=resnrm,"solve_ok"=solve_ok))
}


feval <- function(file.name, ...){
  return(do.call(file.name,list(...)))
}

precondfun <- function(par,r){
  precond <- 0
  if ("precond" %in% names(par)){precoind <- par$precond}
  if (precond == 0){
    q <- r
  }
  else if (precond == 1){
    q <- par$invdiagM * r
  }
  else if (precond == 2){
    dlast  <- par$d[length(par$d)]
    temp1 <- 2/dlast
    temp2 <- 1/par$d - temp1
    q <- temp1 * r + par$V %*% (temp2 * (par$Vt*r))
  }
  else if (precond == 3){
    q <- par$inv
  }
  else if (precond == 4){
    q <- backsolve(par$L,r)
  }
  
  return(q)
}
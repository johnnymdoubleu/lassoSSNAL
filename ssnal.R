library(readr)
library(rmatio)
library(Rcpp)
library(RSpectra)
library(Matrix)

source("folder/Classic_Lasso_SSNAL.R")
source("folder/Classic_Lasso_SSNAL_main.R")
source("folder/Classic_Lasso_SSNCG.R")
source("folder/proj_inf.R")
source("folder/linsyssolve.R")
source("folder/findstep.R")
source("folder/psqmry.R")
source("folder/matvec_ClassicLasso.R")
source("folder/findnnz.R")
sourceCpp("folder/mex_matrix_mult.cpp")
sourceCpp("folder/mexsigma_update_classic_Lasso_SSNAL.cpp")

lipfun <- function(b, A){
  return(t(t(A%*%b) %*% A))
}


#A : a design matrix containing explanatory variable.
#b : a vector containing response variable.
#lambda : a vector containing range of lambda in the power of 10.
#        eg) For warmstarting please input as a vector of length 2.
#            c(0, -5), where it will turn to 10^seq(0, -5, length = 10)
#            Without warmstarting please input as a numerical variable
#            (Note that the larger power must come first in the vector)
#stoptol : a numerical variable corresponding to the stopping tolerance
#          The default value is 1e-6.
#maxiter : a numerical variable corresponding to the maximum number of 
#         iteration. The default value is 1000.
#warmstart : a boolean variable corresponding to perform warmstarting.
#            The default value is FALSE.
ssnal <- function(A, b, lambda = c(0, -5), stoptol = 1e-6, printyes=TRUE,
                  maxiter = 1000, warmstart = FALSE){

  eps <- 2.220446e-16 # Copy the MATLAB eps essentially 
  n <- ncol(A)
  c <- lambda[1]
  # c <- 10^(-3)
  
  maxA <- max(abs(t(t(b) %*% A)))
  eigs_AtA <- eigs_sym(lipfun, k = 1, n = n, args = A)
  Lip <- eigs_AtA$values
  stoptol <- 1e-6
  opts <- list(
    stoptol = stoptol,
    Lip = Lip,
    Ascale = 1,
    maxiter = maxiter,
    printyes = printyes
  )
  
  if (warmstart) {
    output.list <- matrix(nrow=9, ncol=10)
    grid <- 10^seq(lambda[1], lambda[2], length = 10)
    x0 <- NULL
    for(i in grid) {
      cat("SSNAL with lambda:", i, "\n")
      rho <- i * maxA
      Rprof()
      clo <- Classic_Lasso_SSNAL(A, b, n, rho, opts, x = x0)
      Rprof(NULL)
      x0 <- clo$info$x
      print("-------------------------")
      warmstart.list <- list(
        lambda = i,
        rho = rho,
        total.iter = clo$info$iter,
        primal.obj = clo$obj[1],
        dual.obj = clo$obj[2],
        min.x = clo$info$min.x,
        max.x = clo$info$max.x,
        nnz = findnnz(clo$info$x,0.999)$k,
        time.taken = summaryRprof()$sampling.time
      )
      # print(do.call(rbind, warmstart.list))
      output.list[,which(grid==i)] <- do.call(rbind, warmstart.list)
    }
    
    row.names(output.list) <- row.names(do.call(rbind, warmstart.list))
    options(scipen = 10)
    print(output.list)
    options(scipen = 0)
    minlambda <- grid[which(output.list[4,]==min(output.list[4,]))]
    cat("Lowest Objective Value = ", min(output.list[4,]), "\n", 
        "lambda value = ", minlambda, "\n")
    cat("-----------------------------------------", "\n")
    cat("Starting with an optimised lambda value")
    # clo <- Classic_Lasso_SSNAL(A, b, n, (minlambda * maxA), opts, x = x0)
    # return(output.list)
  }
  else{
    # Rprof()
    rho <- c * maxA
    clo <- Classic_Lasso_SSNAL(A, b, n, rho, opts)
    # Rprof(NULL)
  }
}



library(readr)
library(rmatio)
library(Rcpp)
library(RSpectra)
library(Matrix)
library(glmnet)
library(tidyverse)

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



ssnal <- function(A, b, lambda, stoptol = 1e-6, maxiter = 1000,
                  warmstart = FALSE){

  eps <- 2.220446e-16 # Copy the MATLAB eps essentially 
  n <- ncol(A)
  c <- lambda
  # c <- 10^(-3)
  
  maxA <- max(abs(t(t(b) %*% A)))
  eigs_AtA <- eigs_sym(lipfun, k = 1, n = n, args = A)
  Lip <- eigs_AtA$values
  stoptol <- 1e-6
  opts <- list(
    stoptol = stoptol,
    Lip = Lip,
    Ascale = 1,
    maxiter = maxiter
  )
  # opts$stoptol <- stoptol
  # opts$Lip <- Lip
  # opts$Ascale <- 1
  # opts$maxiter <- 10
  
  if (warmstart) {
    output.list <- matrix(nrow=9, ncol=10)
    grid <- 10^seq(-5, 0, length = 10)
    for(lambda in grid) {
      rho <- lambda * maxA
      Rprof()
      clo <- Classic_Lasso_SSNAL(A, b, n, rho, opts)
      Rprof(NULL)
      print("-------------------------")
      warmstart.list <- list(
        lambda = lambda,
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
      output.list[,which(grid==lambda)] <- do.call(rbind, warmstart.list)
    }
    
    row.names(output.list) <- row.names(do.call(rbind, warmstart.list))
    options(scipen = 10)
    output.list
    options(scipen = 0)
    cat("Lowest Objective Value = ", min(output.list[4,]), "\n", 
        "lambda value = ", grid[which(output.list[4,]==min(output.list[4,]))], "\n")
  }
  else{
    Rprof()
    clo <- Classic_Lasso_SSNAL(A, b, n, rho, opts)
    Rprof(NULL)
  }
}



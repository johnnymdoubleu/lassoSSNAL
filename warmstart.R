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


# data <- read.mat("UCIdata/abalone_scale_expanded7.mat")    #lassoSSNAL
# data <- read.mat("UCIdata/space_ga_scale_expanded9.mat")   #lassoSSNAL
# data <- read.mat("UCIdata/bodyfat_scale_expanded7.mat")    #lassoSSNAL
# data <- read.mat("UCIdata/pyrim_scale_expanded5.mat")      #lassoSSNAL
# data <- read.mat("UCIdata/housing_scale_expanded7.mat")    #lassoSSNAL
# data <- read.mat("UCIdata/triazines_scale_expanded4.mat")  #lassoSSNAL
data <- read.mat("UCIdata/mpg_scale_expanded7.mat")        #lassoSSNAL

lipfun <- function(b, A){
  return(t(t(A%*%b) %*% A))
}

A <- data$A
b <- data$b

eps <- 2.220446e-16 # Copy the MATLAB eps essentially 
n <- ncol(A)
c <- 10^(-3)


maxA <- max(abs(t(t(b) %*% A)))
eigs_AtA <- eigs_sym(lipfun, k = 1, n = n, args = A)
Lip <- eigs_AtA$values
stoptol <- 1e-6
opts <- c()
opts$stoptol <- stoptol
opts$Lip <- Lip
opts$Ascale <- 1
# opts$maxiter <- 10

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
    min.x = clo$info$minx,
    max.x = clo$info$max,
    nnz = findnnz(clo$info$x,0.999)$k,
    time.taken = summaryRprof()$sampling.time
  )
  output.list[,which(grid==lambda)] <- do.call(rbind, warmstart.list)
}

row.names(output.list) <- row.names(do.call(rbind, warmstart.list))
options(scipen = 10)
output.list
options(scipen = 0)
cat("Lowest Objective Value = ", min(output.list[4,]), "\n", 
    "lambda value = ", grid[which(output.list[4,]==min(output.list[4,]))], "\n")

# cat("Primal Objective Value = ", clo$obj[1], "\n")
# cat("Dual Objective Value = ", clo$obj[2], "\n")
# cat("min(X) = ", clo$info$minx, "\n")
# cat("max(X) = ", clo$info$max, "\n")
# cat("nnz = ", findnnz(clo$info$x,0.999)$k, "\n")



########## Only for Methylation ##############
# A <- read_delim("UCIdata/GSE40279_average_beta.txt", "\t", col_names = TRUE)
# A[,1] <- NULL
# A <- as.matrix(A)
# b <- as.vector(read.csv("UCIdata/sample.csv", header=FALSE)[,3])
# A <- t(A)              
##############################################

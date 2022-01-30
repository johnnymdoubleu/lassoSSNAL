library(Matrix)
# library(R.matlab)
library(rmatio)
library(RSpectra)
library(Rcpp)

source("lassoSSNAL/Classic_Lasso_SSNAL.R")
source("lassoSSNAL/Classic_Lasso_SSNAL_main.R")
source("lassoSSNAL/Classic_Lasso_SSNCG.R")
source("lassoSSNAL/proj_inf.R")

sourceCpp("lassoSSNAL/test.cpp")
# microbenchmark( (t(data$A)%*%data$A), eigenMatMult(t(data$A),data$A, n_cores=4))

eps <- 2.220446e-16 # Copy the MATLAB eps essentially

#30secs
Rprof()
# data <- readMat("UCIdata/abalone_scale_expanded7.mat")
data <- read.mat("UCIdata/abalone_scale_expanded7.mat")

# AtA <- t(data$A) %*% data$A #2mins15secs

# AtA <- eigenMatMult(t(data$A),data$A, n_cores=4)
AtA <- eigenTransMatMult(data$A, n_cores=4)
eigs_AtA <- eigs_sym(AtA,1) #instant


#30secs

A <- data$A
At <- t(data$A)
b <- data$b

c <- 10^(-4) ## THIS IS LAMBDA
rho <- c*max(abs(At%*%b))
# rho <- c*max(abs(eigenMatMult(At, b, n_cores=4)))
Lip <- eigs_AtA$values
stoptol <- 1e-6
n <- ncol(A)

opts <- c()
opts$stoptol <- stoptol
opts$Lip <- Lip
opts$Ascale <- 1
# Rprof(NULL)
# summaryRprof()
#2mins
test<-Classic_Lasso_SSNAL(A,b,n,rho,opts)
test
Rprof(NULL)
summaryRprof()

#as.numeric(strsplit(format(Sys.time(), "%Y %m %d %H %M %S")," ")[[1]])/rep(1000,6)
#A%*%Diagonal(x=test$dscale)

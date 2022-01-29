library(Matrix)
library(R.matlab)
library(RSpectra)
<<<<<<< Updated upstream
=======
library(Matrix)
library(Rcpp)
library(microbenchmark)
>>>>>>> Stashed changes

source("lassoSSNAL/Classic_Lasso_SSNAL.R")
source("lassoSSNAL/Classic_Lasso_SSNAL_main.R")
source("lassoSSNAL/Classic_Lasso_SSNCG.R")
source("lassoSSNAL/proj_inf.R")
<<<<<<< Updated upstream

eps <- 2.220446e-16 # Copy the MATLAB eps essentially

data <- readMat("UCIdata/abalone_scale_expanded7.mat")
Rprof()
AtA <- t(data$A) %*% data$A #2mins15secs
summaryRprof()
=======

sourceCpp("lassoSSNAL/test.cpp")
# microbenchmark( (t(data$A)%*%data$A), eigenMatMult(t(data$A),data$A, n_cores=4))

eps <- 2.220446e-16 # Copy the MATLAB eps essentially

#30secs
Rprof()
data <- readMat("UCIdata/abalone_scale_expanded7.mat")
# AtA <- t(data$A) %*% data$A #2mins15secs

AtA <- eigenMatMult(t(data$A),data$A, n_cores=4)

>>>>>>> Stashed changes
eigs_AtA <- eigs_sym(AtA,1) #instant
Rprof(NULL)

#30secs

A <- data$A
At <- t(data$A)
b <- data$b

c <- 10^(-4) ## THIS IS LAMBDA
rho <- c*max(abs(eigenMatMult(At, b, n_cores=4)))
Lip <- eigs_AtA$values
stoptol <- 1e-6
n <- ncol(A)

opts <- c()
opts$stoptol <- stoptol
opts$Lip <- Lip
opts$Ascale <- 1

#2mins
test<-Classic_Lasso_SSNAL(A,b,n,rho,opts)
Rprof(NULL)
summaryRprof()
#as.numeric(strsplit(format(Sys.time(), "%Y %m %d %H %M %S")," ")[[1]])/rep(1000,6)
#A%*%Diagonal(x=test$dscale)

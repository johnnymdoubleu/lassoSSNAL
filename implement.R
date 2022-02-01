#library(R.matlab)
library(rmatio)
library(Rcpp)
library(RSpectra)
library(Matrix)
library(DWDLargeR)

source("Classic_Lasso_SSNAL.R")
source("Classic_Lasso_SSNAL_main.R")
source("Classic_Lasso_SSNCG.R")
source("proj_inf.R")

sourceCpp("test.cpp")
sourceCpp("mexsigma_update_classic_Lasso_SSNAL.cpp")

eps <- 2.220446e-16 # Copy the MATLAB eps essentially

#30secs
# data <- readMat("UCIdata/abalone_scale_expanded7.mat")
data <- read.mat("UCIdata/abalone_scale_expanded7.mat")

# AtA <- t(data$A) %*% data$A #2mins15secs

# AtA <- eigenMatMult(t(data$A),data$A, n_cores=4)
AtA <- eigenTransMatMult(data$A, n_cores=4)

eigs_AtA <- eigs_sym(AtA,1) #instant

A <- data$A
At <- t(data$A)
b <- data$b

c <- 10^(-4) ## THIS IS LAMBDA
rho = c*max(abs(At%*%b))
Lip <- eigs_AtA$values
stoptol <- 1e-6
n <- ncol(A)

opts <- c()
opts$stoptol <- stoptol
opts$Lip <- Lip
opts$Ascale <- 1

#2mins
test<-Classic_Lasso_SSNAL(A,b,n,rho,opts)

#as.numeric(strsplit(format(Sys.time(), "%Y %m %d %H %M %S")," ")[[1]])/rep(1000,6)
#A%*%Diagonal(x=test$dscale)


source("Classic_Lasso_SSNAL.R")
source("Classic_Lasso_SSNAL_main.R")
source("Classic_Lasso_SSNCG.R")
source("proj_inf.R")
source("linsyssolve.R")
source("findstep.R")
source("psqmry.R")
source("matvec_ClassicLasso.R")

Classic_Lasso_SSNAL(A,b,n,rho,opts)
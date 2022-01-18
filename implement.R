library(Matrix)
library(R.matlab)
library(RSpectra)

source("lassoSSNAL/Classic_Lasso_SSNAL.R")
source("lassoSSNAL/Classic_Lasso_SSNAL_main.R")
source("lassoSSNAL/Classic_Lasso_SSNCG.R")
source("lassoSSNAL/proj_inf.R")

eps <- 2.220446e-16 # Copy the MATLAB eps essentially

data <- readMat("UCIdata/abalone_scale_expanded7.mat")
Rprof()
AtA <- t(data$A) %*% data$A #2mins15secs
summaryRprof()
eigs_AtA <- eigs_sym(AtA,1) #instant
Rprof(NULL)

#30secs

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

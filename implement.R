#library(R.matlab)
library(rmatio)
library(Rcpp)
library(RSpectra)
library(Matrix)
# library(DWDLargeR)

source("lassoSSNAL/Classic_Lasso_SSNAL.R")
source("lassoSSNAL/Classic_Lasso_SSNAL_main.R")
source("lassoSSNAL/Classic_Lasso_SSNCG.R")
source("lassoSSNAL/proj_inf.R")
source("lassoSSNAL/linsyssolve.R")
source("lassoSSNAL/findstep.R")
source("lassoSSNAL/psqmry.R")
source("lassoSSNAL/matvec_ClassicLasso.R")
source("lassoSSNAL/findnnz.R")
sourceCpp("lassoSSNAL/test.cpp")
sourceCpp("lassoSSNAL/mexsigma_update_classic_Lasso_SSNAL.cpp")

# Rprof()
eps <- 2.220446e-16 # Copy the MATLAB eps essentially
# data <- readMat("UCIdata/abalone_scale_expanded7.mat")
#30secs

# data <- read.mat("UCIdata/abalone_scale_expanded7.mat")    #working
# data <- read.mat("UCIdata/space_ga_scale_expanded9.mat")   #working
# data <- read.mat("UCIdata/bodyfat_scale_expanded7.mat")
data <- read.mat("UCIdata/pyrim_scale_expanded5.mat")
# data <- read.mat("UCIdata/housing_scale_expanded7.mat")
# data <- read.mat("UCIdata/triazines_scale_expanded4.mat")
# data <- read.mat("UCIdata/mpg_scale_expanded7.mat")        #working


A <- data$A
At <- t(data$A)
AtA <- eigenMapMatMult(At, A, n_cores=4)
eigs_AtA <- eigs_sym(AtA,1) #instant

b <- data$b

c <- 10^(-4) ## THIS IS LAMBDA
# rho = c*max(abs(At%*%b))
rho <- c*max(abs(eigenVecMatMult(At, b,4)))
Lip <- eigs_AtA$values
stoptol <- 1e-6
n <- ncol(A)

opts <- c()
opts$stoptol <- stoptol
opts$Lip <- Lip
opts$Ascale <- 1

#as.numeric(strsplit(format(Sys.time(), "%Y %m %d %H %M %S")," ")[[1]])/rep(1000,6)
Rprof()
# clo <- Classic_Lasso_SSNAL(A,b,n,rho,opts)
Rprof(NULL)
summaryRprof()
print("-------------------------")
cat("min(X) = ",clo$info$minx,"\n")
cat("max(X) = ",clo$info$max,"\n")
cat("nnz = ",findnnz(clo$info$x,0.999)$k,"\n")

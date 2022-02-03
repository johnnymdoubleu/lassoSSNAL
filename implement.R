library(Matrix)
# library(R.matlab)
library(rmatio)
library(RSpectra)
library(Rcpp)

source("lassoSSNAL/Classic_Lasso_SSNAL.R")
source("lassoSSNAL/Classic_Lasso_SSNAL_main.R")
source("lassoSSNAL/Classic_Lasso_SSNCG.R")
source("lassoSSNAL/linsyssolve.R")
source("lassoSSNAL/findstep.R")
source("lassoSSNAL/proj_inf.R")
source("lassoSSNAL/psqmry.R")
source("lassoSSNAL/matvec_ClassicLasso.R")
source("lassoSSNAL/findnnz.R")

sourceCpp("lassoSSNAL/test.cpp")
sourceCpp("lassoSSNAL/mexsigma_update_classic_Lasso_SSNAL.cpp")

eps <- 2.220446e-16 # Copy the MATLAB eps essentially

#30secs
Rprof()
# data <- readMat("UCIdata/abalone_scale_expanded7.mat")
data <- read.mat("UCIdata/abalone_scale_expanded7.mat")

# newA <- Matrix(data$A, sparse = TRUE)
A <- data$A
At <- t(data$A)
b <- data$b
AtA <- eigenMapMatMult(At, A,4)
# AtA <- eigenTransMapMatMult(data$A, n_cores=8)
# AtA <- crossCpp(data$A)


eigs_AtA <- eigs_sym(AtA,1) #instant


#30secs


# Rprof(NULL)
# summaryRprof()

c <- 10^(-4) ## THIS IS LAMBDA
# rho <- c*max(abs(At%*%b))
rho <- c*max(abs(eigenMapMatMult(At, b, n_cores=4)))
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
# Rprof()
clo <- Classic_Lasso_SSNAL(A,b,n,rho,opts)

Rprof(NULL)
summaryRprof()

# library(profr)
# ggplot.profr(parse_rprof("profile1.out"))

print("-------------------------")
cat("min(X) = ",clo$info$minx,"\n")
cat("max(X) = ",clo$info$max,"\n")
cat("nnz = ",findnnz(clo$info$x,0.999)$k,"\n")
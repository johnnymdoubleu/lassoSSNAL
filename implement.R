# library(data.table)
library(readr)
library(rmatio)
library(Rcpp)
library(RSpectra)
library(Matrix)
library(glmnet)
# library(DWDLargeR)

source("Classic_Lasso_SSNAL.R")
source("Classic_Lasso_SSNAL_main.R")
source("Classic_Lasso_SSNCG.R")
source("proj_inf.R")
source("linsyssolve.R")
source("findstep.R")
source("psqmry.R")
source("matvec_ClassicLasso.R")
source("findnnz.R")
sourceCpp("test.cpp")
sourceCpp("mexsigma_update_classic_Lasso_SSNAL.cpp")

# Rprof()
eps <- 2.220446e-16 # Copy the MATLAB eps essentially
#30secs

# data <- read.mat("UCIdata/abalone_scale_expanded7.mat")    #working
# data <- read.mat("UCIdata/space_ga_scale_expanded9.mat")   #working
# data <- read.mat("UCIdata/bodyfat_scale_expanded7.mat")    #working
# data <- read.mat("UCIdata/pyrim_scale_expanded5.mat")      #working
 data <- read.mat("UCIdata/housing_scale_expanded7.mat")    #working
# data <- read.mat("UCIdata/triazines_scale_expanded4.mat")  #working
# data <- read.mat("UCIdata/mpg_scale_expanded7.mat")        #working

# data <- read.mat("UCIdata/E2006.train.mat")
# data <- read.mat("UCIdata/E2006.test.mat")

 
 
 
A <- read_delim("GSE40279_average_beta.txt", "\t", col_names = TRUE)
# A <- A[-1,]
A[,1] <- NULL
A <- as.matrix(A)
b <- as.vector(read.csv("sample.csv", header=FALSE)[,3])




lipfun <- function(b, A){
  return(t(t(A%*%b) %*% A))
}


#A <- t(A)

# A <- as.matrix(data$A)
# b <- data$b

 A <- data$A
 b <- data$b
 
 
 
grid <- 10^seq(-6, -20, length = 10) * max(abs(t(t(b) %*% A)))
glambda <- cv.glmnet(A,b,lambda = grid)
# glmnet(A,b,alpha=1, lambda = grid)
plot(glambda)
 
 
n <- ncol(A)

c <- 10^(-3) ## THIS IS LAMBDA
#c <- 3.727594e-03
#c <- glambda$lambda.min
rho <- c * max(abs(t(t(b) %*% A)))
rho <- 0.04186
# Rprof(NULL)
# summaryRprof()

eigs_AtA <- eigs_sym(lipfun, k = 1, n = n, args = A)
Lip <- eigs_AtA$values
stoptol <- 1e-1


opts <- c()
opts$stoptol <- stoptol
opts$Lip <- Lip
opts$Ascale <- 1

#as.numeric(strsplit(format(Sys.time(), "%Y %m %d %H %M %S")," ")[[1]])/rep(1000,6)
Rprof()
clo <- Classic_Lasso_SSNAL(A, b, n, rho, opts)
Rprof(NULL)
summaryRprof()
print("-------------------------")
cat("min(X) = ", clo$info$minx, "\n")
cat("max(X) = ", clo$info$max, "\n")
cat("nnz = ", findnnz(clo$info$x,0.999)$k, "\n")


coefs <- clo$x
objective <- 0.5*sum((b - (A %*% coefs))^2) + 
  rho * sum(abs(coefs))
objective


#c<-1.14016
Rprof()
test <- glmnet(A,b/sd(b),alpha=1,lambda=rho/(sd(b)*length(b)),intercept=FALSE,standardize=FALSE)
Rprof(NULL)
summaryRprof()

coefs <- coef(test)[-1]*sd(b)
objective <- 0.5*sum((b - (A %*% coefs))^2) + 
  rho * sum(abs(coefs))
objective
min(coef(test)[-1]*sd(b))
max(coef(test)[-1]*sd(b))
findnnz(as.matrix(coef(test)[-1]*sd(b)),0.999)$k



plot(A%*%coef(test)[-1],b)
plot(predict(test,A),b)
plot(A%*%clo$x,b)

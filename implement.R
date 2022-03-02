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
 
 
 
grid <- 10^seq(-6, -5, length = 4) * max(abs(t(t(b) %*% A)))
glambda <- cv.glmnet(A,b,lambda = grid)
# glmnet(A,b,alpha=1, lambda = grid)
plot(glambda)
 
 
n <- ncol(A)

c <- 10^(-4) ## THIS IS LAMBDA
#c <- 3.727594e-03
#c <- glambda$lambda.min
rho <- c * max(abs(t(t(b) %*% A)))





rho <- 0.04186*100
# Rprof(NULL)
# summaryRprof()

eigs_AtA <- eigs_sym(lipfun, k = 1, n = n, args = A)



Lip <- eigs_AtA$values
stoptol <- 0.99
opts <- c()
opts$stoptol <- stoptol
opts$Lip <- Lip
opts$Ascale <- 1
opts$maxiter <- 10

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















##### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##### CROSS-VALIDATION CODE FOR METHYLATION DATA
##### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





# Make consistent groups
set.seed(1)
numobs<-656
nx<-sample(1:numobs)
levels<-10
gmem<-split(nx,nx%%levels)
#groups[1]<-nx[as.numeric(unlist(gmem[1]))]
#groups <- c()
#for(i in 1:10) {
#  groups <- c(groups,nx[as.numeric(unlist(gmem[i]))])
#}




rhomult <- 31.62
rho <- 0.04186*rhomult

mses <- c()
objs <- c()
nnzs <- c()


Lips <- readRDS("Lips.rds")

#test_gi <- 2
for(test_gi in 1:10) {
  cat("\n\nNow doing group #",test_gi,"\n\n")
  
  tg_mems <- nx[as.numeric(unlist(gmem[test_gi]))]
  test_A <- A[tg_mems,]
  test_b <- b[tg_mems]
  
  train_A <- A[-tg_mems,]
  train_b <- b[-tg_mems]
  
  
  #eigs_AtA <- eigs_sym(lipfun, k = 1, n = n, args = train_A)
  #Lip <- eigs_AtA$values
  Lip <- Lips[test_gi]
  stoptol <- 0.99
  opts <- c()
  opts$stoptol <- stoptol
  opts$Lip <- Lip
  opts$Ascale <- 1
  opts$maxiter <- 10
  
  clo <- Classic_Lasso_SSNAL(train_A, train_b, n, rho, opts)
  
  #plot(train_A%*%clo$x,train_b)
  mse <- sum((test_b - test_A%*%clo$x)^2)/length(test_b)
  nnz <- findnnz(clo$info$x,0.999)$k

  coefs <- clo$x
  obj <- 0.5*sum((train_b - (train_A %*% coefs))^2) + 
    rho * sum(abs(coefs))
  
  
  
  ### Append vars to vectors this rhomult
  
  mses <- c(mses, mse)
  objs <- c(objs, obj)
  nnzs <- c(nnzs, nnz)
}



write.csv(nnzs,paste0(gsub('\\.','-',toString(rhomult)),"_nnz.csv"))
write.csv(objs,paste0(gsub('\\.','-',toString(rhomult)),"_obj.csv"))
write.csv(mses,paste0(gsub('\\.','-',toString(rhomult)),"_mse.csv"))
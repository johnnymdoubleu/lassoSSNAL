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
 data <- read.mat("UCIdata/space_ga_scale_expanded9.mat")   #working
# data <- read.mat("UCIdata/bodyfat_scale_expanded7.mat")    #working
# data <- read.mat("UCIdata/pyrim_scale_expanded5.mat")      #working
# data <- read.mat("UCIdata/housing_scale_expanded7.mat")    #working
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


A <- t(A)

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
stoptol <- 1e-6
opts <- c()
opts$stoptol <- stoptol
opts$Lip <- Lip
opts$Ascale <- 1
opts$maxiter <- 5000

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




rhomult <- 100
rho <- 0.04186*rhomult

mses <- c()
objs <- c()
nnzs <- c()


Lips <- readRDS("Lips.rds")

#test_gi <- 2
for(test_gi in 8:10) {
  cat("\n\nNow doing group #",test_gi,"\n\n")
  
  tg_mems <- nx[as.numeric(unlist(gmem[test_gi]))]
  test_A <- A[tg_mems,]
  test_b <- b[tg_mems]
  
  train_A <- A[-tg_mems,]
  train_b <- b[-tg_mems]
  
  
  #eigs_AtA <- eigs_sym(lipfun, k = 1, n = n, args = train_A)
  #Lip <- eigs_AtA$values
  Lip <- Lips[test_gi]
  stoptol <- 0.01
  opts <- c()
  opts$stoptol <- stoptol
  opts$Lip <- Lip
  opts$Ascale <- 1
  opts$maxiter <- 20
  
  clo <- Classic_Lasso_SSNAL(train_A, train_b, n, rho, opts)
  
  #plot(train_A%*%clo$x,train_b)
  mse <- sum((test_b - test_A%*%clo$x)^2)/length(test_b)
  nnz <- findnnz(clo$info$x,0.999)$k

  coefs <- clo$x
  obj <- 0.5*sum((train_b - (train_A %*% coefs))^2) + 
    rho * sum(abs(coefs))
  
  
  cat("\nmse=",mse,"\n")
  
  ### Append vars to vectors this rhomult
  
  mses <- c(mses, mse)
  objs <- c(objs, obj)
  nnzs <- c(nnzs, nnz)
}



write.csv(nnzs,paste0(gsub('\\.','-',toString(rhomult)),"_nnz3.csv"))
write.csv(objs,paste0(gsub('\\.','-',toString(rhomult)),"_obj3.csv"))
write.csv(mses,paste0(gsub('\\.','-',toString(rhomult)),"_mse3.csv"))




















#glfits <- cv.glmnet(A,b,alpha=1,lambda=rhos_cvg*10/(9*length(b)),standardize=FALSE,intercept=FALSE)









r_ssnal_cv_df <- data.frame(log_lam=numeric(),mse=numeric(),ciw=numeric())
r_ssnal_mults <- c(100,316.2,1000,3162,10000)
for(tmi in r_ssnal_mults) {
  tll <- log(tmi*0.04816*9/(10*656))
  tmi_fsf <- paste0(gsub('\\.','-',toString(tmi)),"_mse.csv")
  tmi_rmd <- read.csv(tmi_fsf)
  mse <- mean(tmi_rmd$x)
  r_ssnal_cv_df <- rbind(r_ssnal_cv_df,c(log_lam=tll,mse=mse,ciw=sd(tmi_rmd$x)/sqrt(10)))
}
colnames(r_ssnal_cv_df)<-c("log_lam","mse","ciw")

glmnet_cv_obj <- readRDS('glfits.rds')
glm_cv_df <- data.frame(log_lam=log(glmnet_cv_obj$lambda),mse=glmnet_cv_obj$cvm)

plot(glm_cv_df$log_lam,glm_cv_df$mse,col='red',pch=19,ylim=c(40,170))
points(r_ssnal_cv_df$log_lam,r_ssnal_cv_df$mse,col='blue',pch=19)

ndp<-7
ggpdf <- data.frame(ll=glm_cv_df$log_lam,mse_ssnal=numeric(ndp),ciw_ssnal=numeric(ndp),mse_glm=numeric(ndp),ciw_glm=numeric(ndp))

ggpdf[1:5,"mse_ssnal"] <- rev(r_ssnal_cv_df$mse)
ggpdf[1:5,"ciw_ssnal"] <- rev(r_ssnal_cv_df$ciw)
ggpdf$mse_glm <- glmnet_cv_obj$cvm
ggpdf$ciw_glm <- glmnet_cv_obj$cvsd
ggpdf[ggpdf==0] <- NA

ggplot(ggpdf,aes(x=ll,y=mse_ssnal))+
  geom_point(aes(x=ll,y=mse_ssnal),size=2.2,col='blue')+
  geom_point(aes(x=ll,y=mse_glm),size=2.2,col='red')+theme_light()

ggpdf2<-data.frame(ll=numeric(12),mse=numeric(12),ciw=numeric(12),grp=factor(12),cl=factor(12))
ggpdf2[1:7,'ll'] <- glm_cv_df$log_lam
ggpdf2[8:12,'ll'] <- glm_cv_df$log_lam[1:5]

ggpdf2[1:7,'mse'] <- ggpdf[1:7,'mse_glm']
ggpdf2[8:12,'mse'] <- ggpdf[1:5,'mse_ssnal']

ggpdf2[1:7,'ciw'] <- ggpdf[1:7,'ciw_glm']
ggpdf2[8:12,'ciw'] <- ggpdf[1:5,'ciw_ssnal']

ggpdf2$grp <- c(rep('glmnet',7),rep('ssnal',5))
ggpdf2$cl <- c(rep('red',7),rep('blue',5))

ggplot(ggpdf2,aes(x=ll,y=mse,group=grp,color=grp))+
  geom_point(size=2.2)+
  geom_errorbar(aes(ymin=mse-ciw, ymax=mse+ciw), width=.2,
                position=position_dodge(0.05))+
  theme_light()+scale_color_manual(values=c("blue","red"))











pdf <- readRDS("cvoplam_predac.rds")
plot(pdf$predicted,pdf$actual,
     xlab="Predicted age",
     ylab="Actual age",cex=0.7)
abline(c(0,1))
















pred_df<-data.frame(predicted=numeric(),actual=numeric())

rhomult <- 1000
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
  stoptol <- 0.01
  opts <- c()
  opts$stoptol <- stoptol
  opts$Lip <- Lip
  opts$Ascale <- 1
  opts$maxiter <- 20
  
  clo <- Classic_Lasso_SSNAL(train_A, train_b, n, rho, opts)
  
  #plot(train_A%*%clo$x,train_b)
  mse <- sum((test_b - test_A%*%clo$x)^2)/length(test_b)
  nnz <- findnnz(clo$info$x,0.999)$k
  
  coefs <- clo$x
  obj <- 0.5*sum((train_b - (train_A %*% coefs))^2) + 
    rho * sum(abs(coefs))
  
  
  cat("\nmse=",mse,"\n")
  
  ### Append vars to vectors this rhomult
  
  mses <- c(mses, mse)
  objs <- c(objs, obj)
  nnzs <- c(nnzs, nnz)
  
  
  
  
  pred_df_sub <- data.frame(predicted=test_A%*%clo$x,actual=test_b)
  pred_df <- rbind(pred_df,pred_df_sub)
  colnames(pred_df)<-c("predicted","actual")
}

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
c <- 10^(-3)                          # THIS IS Lambda_c
#c <- 3.727594e-03
#c <- glambda$lambda.min
rho <- c * max(abs(t(t(b) %*% A)))    # This is Lambda
# rho <- 0.04186*100

eigs_AtA <- eigs_sym(lipfun, k = 1, n = n, args = A)
Lip <- eigs_AtA$values
# stoptol <- 0.99
stoptol <- 1e-6
opts <- c()
opts$stoptol <- stoptol
opts$Lip <- Lip
opts$Ascale <- 1
# opts$maxiter <- 10

Rprof()
clo <- Classic_Lasso_SSNAL(A, b, n, rho, opts)
Rprof(NULL)
# summaryRprof()
print("-------------------------")
cat("min(X) = ", clo$info$minx, "\n")
cat("max(X) = ", clo$info$max, "\n")
cat("nnz = ", findnnz(clo$info$x,0.999)$k, "\n")




########## Only for Methylation ##############
# A <- read_delim("UCIdata/GSE40279_average_beta.txt", "\t", col_names = TRUE)
# A[,1] <- NULL
# A <- as.matrix(A)
# b <- as.vector(read.csv("UCIdata/sample.csv", header=FALSE)[,3])
# A <- t(A)              
##############################################


# 
# 
# coefs <- clo$x
# objective <- 0.5*sum((b - (A %*% coefs))^2) + 
#   rho * sum(abs(coefs))
# objective
# 
# # c<-1.14016
# Rprof()
# test <- glmnet(A,b/sd(b),alpha=1,lambda=rho/(sd(b)*length(b)),intercept=FALSE,standardize=FALSE)
# Rprof(NULL)
# summaryRprof()
# 
# coefs <- coef(test)[-1]*sd(b)
# objective <- 0.5*sum((b - (A %*% coefs))^2) + 
#   rho * sum(abs(coefs))
# objective
# min(coef(test)[-1]*sd(b))
# max(coef(test)[-1]*sd(b))
# findnnz(as.matrix(coef(test)[-1]*sd(b)),0.999)$k

# 
# plot(A%*%coef(test)[-1],b)
# plot(predict(test,A),b)
# plot(A%*%clo$x,b)

##### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##### CROSS-VALIDATION CODE FOR UCI / Statlib Data
##### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Make consistent groups

####
# set.seed(1)
# 
# 
# nx<-sample(1:nrow(A))
# levels<-2
# gmem<-split(nx,nx%%levels)
# 
# nama <- "triazines4"
# grid <- 10^seq(-14/3, -2, length = 10) * max(abs(t(t(b) %*% A)))
# saveRDS(grid, paste0("lassoSSNAL/UCI/",nama,"/", nama, "_lambda.rds"))
# for (i in grid[2:10]){
#   mses.ssnal <- mses.glm <- objs.ssnal <- objs.glm <- nnzs.ssnal <- nnzs.glm <- c()
#   cat("\n\nNow doing rho value", i,"\n\n\n\n\n\n\n\n")
#   for(test_gi in c(1:2)) {
#     cat("\n\nNow doing group #",test_gi,"\n\n")
#   
#     tg_mems <- nx[as.numeric(unlist(gmem[test_gi]))]
#     test_A <- A[tg_mems,]
#     test_b <- b[tg_mems]
#   
#     train_A <- A[-tg_mems,]
#     train_b <- b[-tg_mems]
#   
#   
#     eigs_AtA <- eigs_sym(lipfun, k = 1, n = ncol(train_A), args = train_A)
#     Lip <- eigs_AtA$values
# 
#     stoptol <- 1e-6
#     opts <- c()
#     opts$stoptol <- stoptol
#     opts$Lip <- Lip
#     opts$Ascale <- 1
#     opts$maxiter <- 500
#     # 
#     rho <- i
#     # 1
#     clo <- Classic_Lasso_SSNAL(train_A, train_b, n, rho, opts)
#     obj <- clo$obj[1]
#     
#     # coefs <- coef(test)[-1]*sd(b)
#     # obj <- 0.5*sum((b - (A %*% coefs))^2) +
#     #   rho * sum(abs(coefs))
#     # #plot(train_A%*%clo$x,train_b)
#     mse <- sum((test_b - test_A%*%clo$x)^2)/length(test_b)
# 
#     nnz <- findnnz(clo$info$x,0.999)$k
# 
#     mses.ssnal <- c(mses.ssnal, mse)
#     nnzs.ssnal <- c(nnzs.ssnal, nnz)
#     objs.ssnal <- c(objs.ssnal, obj)
#     
#     
#     test <- glmnet(train_A, train_b/sd(train_b),
#                    alpha = 1,lambda = rho/(sd(train_b)*length(train_b)),
#                    intercept = FALSE, standardize = FALSE)
#     #
#     coefs <- coef(test)[-1]*sd(b)
#   
#     mse <- sum((test_b - test_A%*%coefs)^2)/length(test_b)
#     
#     obj <- 0.5*sum((train_b - (train_A %*% coefs))^2) +
#       rho * sum(abs(coefs))
#     nnz <- findnnz(as.matrix(coefs),0.999)$k
#     min(coef(test)[-1]*sd(b))
#     max(coef(test)[-1]*sd(b))
#     
#     ### Append vars to vectors this rhomult
#     
#     mses.glm <- c(mses.glm, mse)
#     objs.glm <- c(objs.glm, obj)
#     nnzs.glm <- c(nnzs.glm, nnz)
#   }
#   
#   #choose the choice of format to save
#   saveRDS(mses.ssnal, paste0("lassoSSNAL/UCI/", nama, "/",nama, "_ssnal", which(grid==i), "_mse.rds"))
#   saveRDS(mses.glm, paste0("lassoSSNAL/UCI/", nama, "/",nama, "_glmnet", which(grid==i), "_mse.rds"))
#   saveRDS(nnzs.ssnal, paste0("lassoSSNAL/UCI/", nama, "/",nama, "_ssnal", which(grid==i), "_nnz.rds"))
#   saveRDS(nnzs.glm, paste0("lassoSSNAL/UCI/", nama, "/",nama, "_glmnet", which(grid==i), "_nnz.rds"))
#   saveRDS(objs.ssnal, paste0("lassoSSNAL/UCI/", nama, "/",nama, "_ssnal", which(grid==i), "_obj.rds"))
#   saveRDS(objs.glm, paste0("lassoSSNAL/UCI/", nama, "/",nama, "_glmnet", which(grid==i), "_obj.rds"))
#   
#   
#   # write.csv(mses.ssnal, paste0("lassoSSNAL/", nama,"/",nama, "_ssnal", which(grid==i), "_mse.csv"))
#   # write.csv(mses.glm, paste0("lassoSSNAL/", nama,"/",nama, "_glmnet", which(grid==i), "_mse.csv"))
#   # write.csv(nnzs.ssnal, paste0("lassoSSNAL/", nama,"/",nama, "_ssnal", which(grid==i), "_nnz.csv"))
#   # write.csv(nnzs.glm, paste0("lassoSSNAL/", nama,"/",nama, "_glmnet", which(grid==i), "_nnz.csv"))
#   # write.csv(objs.ssnal, paste0("lassoSSNAL/", nama,"/",nama, "_ssnal", which(grid==i), "_obj.csv"))
#   # write.csv(objs.glm, paste0("lassoSSNAL/", nama,"/",nama, "_glmnet", which(grid==i), "_obj.csv"))
# }


###-----------------------------------------###
#################plotting data#################
###-----------------------------------------###

# for (j in c("abalone7","bodyfat7","housing7","mpg7","pyrim5","space_ga9","triazines4")){
#   nama <- j
#   r_ssnal_cv_df <- data.frame(log_lam=numeric(),mse=numeric(), ciw=numeric(),nnz=numeric())
#   glm_cv_df <- data.frame(log_lam=numeric(),mse=numeric(),ciw=numeric(), nnz=numeric())
#   each.lambda <- readRDS(paste0("lassoSSNAL/UCI/",nama,"/", nama, "_lambda.rds"))
# 
#   for (i in each.lambda) {
#     tll <- log(i)
#     num <- which(each.lambda==i)
#     tmi_rmd <- readRDS(paste0("lassoSSNAL/UCI/",nama,"/", nama, "_ssnal", num, "_mse.rds"))
#     tmi_nnz <- readRDS(paste0("lassoSSNAL/UCI/",nama,"/", nama,  "_ssnal", num, "_nnz.rds"))
# 
#     # tmi_rmd <- read.csv(paste0("lassoSSNAL/", nama,"/",nama, "_ssnal", num, "_mse.csv"))
#     # tmi_nnz <- read.csv(paste0("lassoSSNAL/", nama,"/",nama,  "_ssnal", num, "_nnz.csv"))
# 
#     # tmi_rmd <- readRDS("lassoSSNAL/Methylation/resobjs.rds")
#     # tmi_intres <- readRDS("lassoSSNAL/Methylation/intres.rds")
#     # tmi_goodlam.test <- readRDS("lassoSSNAL/Methylation/goodlam_1000_methtest_forobjonly.rds")
#     # tmi_goodlam <- readRDS("lassoSSNAL/Methylation/goodlam_1000_meth_forobjonly.rds")
#     r_ssnal_cv_df <- rbind(r_ssnal_cv_df, c(log_lam = tll,
#                                             mse = mean(tmi_rmd),
#                                             ciw = sd(tmi_rmd)/sqrt(2),
#                                             nnz = round(mean(tmi_nnz))))
#     colnames(r_ssnal_cv_df) <- c("log_lam","mse", "ciw", "nnz")
#     glmnet_cv_obj <- readRDS(paste0("lassoSSNAL/UCI/",nama,"/", nama, "_glmnet", num, "_mse.rds"))
#     glmnet_nnz <- readRDS(paste0("lassoSSNAL/UCI/",nama,"/", nama, "_glmnet", num, "_nnz.rds"))
#     # glmnet_cv_obj <- read.csv(paste0("lassoSSNAL/", nama,"/",nama, "_glmnet", which(grid==i), "_mse.csv"))
#     # glmnet_nnz <- read.csv(paste0("lassoSSNAL/", nama,"/",nama,"_glmnet", which(grid==i), "_nnz.csv"))
# 
#     glm_cv_df <- rbind(glm_cv_df, c(log_lam = tll,
#                                     mse = mean(glmnet_cv_obj),
#                                     ciw = sd(glmnet_cv_obj)/sqrt(2),
#                                     nnz = round(mean(glmnet_nnz))))
#     colnames(glm_cv_df) <- c("log_lam", "mse", "ciw", "nnz")
#   }
#   ndp <- 10
#   ggpdf <- data.frame(ll = glm_cv_df$log_lam,
#                       mse_ssnal=r_ssnal_cv_df$mse,
#                       ciw_ssnal=r_ssnal_cv_df$ciw,
#                       nnz_ssnal=r_ssnal_cv_df$nnz,
#                       mse_glm=glm_cv_df$mse,
#                       ciw_glm=glm_cv_df$ciw,
#                       nnz_glm=glm_cv_df$nnz)
# 
#   ggpdf["mse_ssnal"] <- r_ssnal_cv_df$mse
#   ggpdf["ciw_ssnal"] <- r_ssnal_cv_df$ciw
#   ggpdf["nnz_ssnal"] <- r_ssnal_cv_df$nnz
#   ggpdf["mse_glm"] <- glm_cv_df$mse
#   ggpdf["ciw_glm"] <- glm_cv_df$ciw
#   ggpdf["nnz_glm"] <- glm_cv_df$nnz
#   ggpdf[ggpdf==0] <- NA
# 
# 
#   # plot(A%*%coef(test)[-1],b)
#   # plot(predict(tmi_goodlam,A),b)
#   # plot(A %*% tmi_goodlam$x, b)
# 
#   # preddata <- data.frame(x=A%*%tmi_goodlam$x, y = b)
#   #
#   # ggplot(preddata, aes(x=x, y = y))+
#   #   geom_point(shape=1, size=2.2)+ geom_smooth(method=lm, se=FALSE, color="grey")
#   #   labs(x="Predicted Age", y="Actual Age") +
#   #   theme_light() + scale_color_manual(values=c("blue", "red")) +
#   #   theme(plot.title = element_text(hjust = 0.5), legend.position="none",
#   #         axis.text = element_text(size= 17),
#   #         axis.title = element_text(size = 17))
#   #
#   # ggsave(filename=paste0("lassoSSNAL/Methylation/predict.png"), plot=last_plot())
# 
# 
#   ggplot(ggpdf,aes(x=ll,y=mse_ssnal))+
#     geom_point(aes(x=ll,y=mse_ssnal),size=2.2,col='blue')+
#     geom_point(aes(x=ll,y=mse_glm),size=2.2,col='red')+theme_light()
# 
#   # tmi_rmd <- readRDS("lassoSSNAL/Methylation/resobjs.rds")
#   # tmi_rmd <- arrange(tmi_rmd, mult)
# 
#   # ggpdf2<-data.frame(ll=numeric(20),mse=numeric(20),grp=factor(20),Algorithm=factor(20))
#   # ggpdf2[1:10,'ll'] <- log(tmi_rmd["mult"]*0.04186) - log(656)
#   # ggpdf2[11:20,'ll'] <- log(tmi_rmd["mult"]*0.04186) - log(656)
#   # ggpdf2[1:10,'mse'] <- tmi_rmd["glmnet"]
#   # ggpdf2[11:20,'mse'] <- tmi_rmd["ssnal"]
# 
#   ggpdf2<-data.frame(ll=numeric(20),mse=numeric(20),ciw=numeric(20),nnz=numeric(20),grp=factor(20),Algorithm=factor(20))
#   ggpdf2[1:10,'ll'] <- glm_cv_df$log_lam
#   ggpdf2[11:20,'ll'] <- glm_cv_df$log_lam
#   ggpdf2[1:10,'mse'] <- ggpdf[,'mse_glm']
#   ggpdf2[11:20,'mse'] <- ggpdf[,'mse_ssnal']
#   ggpdf2[1:10,'ciw'] <- ggpdf[,'ciw_glm']
#   ggpdf2[11:20,'ciw'] <- ggpdf[,'ciw_ssnal']
#   ggpdf2[1:10,'nnz'] <- ggpdf[,'nnz_glm']
#   ggpdf2[11:20,'nnz'] <- ggpdf[,'nnz_ssnal']
#   ggpdf2[,"nnz"] <- round(ggpdf2[,"nnz"])
#   ggpdf2$grp <- c(rep('glmnet',10),rep('SSNAL',10))
#   ggpdf2$Algorithm <- c(rep('glmnet',10),rep('SSNAL',10))
#   
#   ggplot(ggpdf2,aes(x = ll,y = mse, group = grp, color = Algorithm))+
#     geom_point(size=2.2) +
#     geom_errorbar(aes(ymin = mse-ciw, ymax = mse+ciw), width = .2)+
#     geom_text(aes(label = nnz), position = position_dodge(width = 1),
#               check_overlap = FALSE, size = 4) +
#     # ggtitle("Objective Values") +
#     # labs(x="log(\u03bb)", y="Objective") +
#     labs(x="log(\u03bb)", y="CV Mean-squared Error") +
#     theme_light() + scale_color_manual(values=c("blue", "red")) +
#     theme(plot.title = element_text(hjust = 0.5, size = 20), legend.position="right",
#           axis.text = element_text(size= 12),
#           axis.title = element_text(size = 20))
# 
#   # ggsave(filename=paste0("lassoSSNAL/Methylation/objval.png"), plot=last_plot())
#   ggsave(filename=paste0("lassoSSNAL/UCI/", nama,".png"), plot=last_plot())
# }

source("lassoSSNAL/ssnal.R", echo = FALSE)

########## UCI Datasets ##############
# data <- read.mat("D:/UCIdata/abalone_scale_expanded7.mat")    #lassoSSNAL
# data <- read.mat("D:/UCIdata/space_ga_scale_expanded9.mat")   #lassoSSNAL
# data <- read.mat("D:/UCIdata/bodyfat_scale_expanded7.mat")    #lassoSSNAL
# data <- read.mat("D:/UCIdata/pyrim_scale_expanded5.mat")      #lassoSSNAL
# data <- read.mat("D:/UCIdata/housing_scale_expanded7.mat")    #lassoSSNAL
# data <- read.mat("D:/UCIdata/triazines_scale_expanded4.mat")  #lassoSSNAL
data <- read.mat("D:/UCIdata/mpg_scale_expanded7.mat")        #lassoSSNAL
A <- data$A
b <- data$b

# Rprof()
# ccc <- ssnal(A,b, lambda = 1e-05, stoptol = 1e-6) #without warmstarting
# Rprof(NULL)
# summaryRprof()
ccc <- ssnal(A,b, lambda = c(0, -5), warmstart = TRUE) #with warmstarting
# ccc <- ssnal(A,b, lambda = 10^(-3), stoptol = 1e-6, printyes = TRUE) #without warmstarting
# ccc <- ssnal(A,b, lambda = c(0, 10^(-3)), warmstart = TRUE, printyes = TRUE, stoptol = 1e-6) #without warmstarting

#---------------------------------------------------------------------
########## Only for Methylation ##############
# A <- read_delim("D:/UCIdata/GSE40279_average_beta.txt", "\t", col_names = TRUE)
# A[,1] <- NULL
# A <- as.matrix(A)
# b <- as.vector(read.csv("UCIdata/sample.csv", header=FALSE)[,3])
# A <- t(A)

# ssnal(A,b, 10^(-3), warmstart = TRUE, maxiter = 15) #with warmstarting
# ssnal(A,b, 10^(-3), warmstart = FALSE, maxiter = 15) #with warmstarting
##############################################
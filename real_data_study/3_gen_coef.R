rm(list=ls())


source("PFAdecomp.R")


library(HiClimR)
library(pracma)



##################################
##################################


data_asian=as.matrix(read.delim("./datafile/data_hc_06.txt", row.names = 1, header= TRUE))

data_asian=as.data.frame(data_asian)

numgene=dim(data_asian)[1]-1

othergene=t(data_asian[2:(numgene+1),])

CHRNA6=t(data_asian[1,])

X <- as.matrix(othergene)
Y <- as.matrix(CHRNA6)
if (nrow(X) != nrow(Y) | ncol(Y) != 1) 
  stop("Dimensions do not match.")

n <- nrow(X)
p <- ncol(X)
Z <- rep(0, p)
for (j in 1:p) {
  #Z[j] <- lm(Y ~ X[, j])$coef[2]
  Z[j] <- lsfit(x = X[, j], y = Y)$coef[2]
}


sigma.est <- sd(Y)
#sd(Y)
#0.331



cormat = fastCor(X,nSplit = 5, upperTri = FALSE, optBLAS = TRUE, verbose = TRUE)
varvec= apply(X, 2, var)

Sigma <- matrix(rep(0, p * p), nrow = p)

for (i in 1:p) {
  for (j in i:p) {
    Sigma[i, j] <-Sigma[j,i]<- ((sigma.est)^2*cormat[i,j])/(n*sqrt(varvec[i])*sqrt(varvec[j]))
  }
}

e=0.01
smallestf=0.9

stdized='Y'

results_PFA=PFAdecomp(Z, Sigma, smallestf,e,stdized)


beta=results_PFA[[1]]

var_K=results_PFA[[2]]

SD=results_PFA[[3]]


path_name="datafile/"

print(path_name)

file_name=paste(path_name,"beta_PFA.csv",sep="/")
write.table(beta, file = file_name, sep = ",", row.names = FALSE,col.names = FALSE )

file_name=paste(path_name,"var_K_PFA.csv",sep="/")
write.table(var_K, file = file_name, sep = ",", row.names = FALSE,col.names = FALSE )

file_name=paste(path_name,"z.csv",sep="/")
write.table(Z, file = file_name, sep = ",", row.names = FALSE,col.names = FALSE )

file_name=paste(path_name,"Sigma_SD.csv",sep="/")
write.table(SD, file = file_name, sep = ",", row.names = FALSE,col.names = FALSE )



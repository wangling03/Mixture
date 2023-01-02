rm(list=ls())

#setwd("E:\\C drive\\Ling MSU paper\\paper 3\\Paper 3\\ZLiao\\real_data_study4")

source("PFAdecomp2.R")

#library(quantreg)
#library(qut)
library(HiClimR)
library(pracma)

#library(varEst)
#remove.packages("varEst")


##################################
##################################



#data_asian=as.matrix(read.delim("E:\\C drive\\Ling MSU paper\\paper 3\\Paper 3\\ZLiao\\real_data_study4\\data_file\\data_hc_05.txt", row.names = 1, header= TRUE))
data_asian=as.matrix(read.delim("./data_file/data_hc_05.txt", row.names = 1, header= TRUE))

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



#sigma.est <- rcv(X, Y)
#0.26

#rcv.est <-sigmarcv(Y, X, cv = FALSE, fit = NA, intercept = TRUE)
#sigma.est=rcv.est$sigmahat
#0.285

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

results_PFA=PFAdecomp2(Z, Sigma, smallestf,e) 


beta=results_PFA[[1]]

var_K=results_PFA[[2]]


path_name="/mnt/home/wangli35/paper_3/real_data_study4/data_file/0.5"

print(path_name)

file_name=paste(path_name,"beta_PFA.csv",sep="/")
write.table(beta, file = file_name, sep = ",", row.names = FALSE,col.names = FALSE )

file_name=paste(path_name,"var_K_PFA.csv",sep="/")
write.table(var_K, file = file_name, sep = ",", row.names = FALSE,col.names = FALSE )

file_name=paste(path_name,"z.csv",sep="/")
write.table(Z, file = file_name, sep = ",", row.names = FALSE,col.names = FALSE )




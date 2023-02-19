rm(list=ls(all=TRUE))



source("GenSigmaMatrix.R")
source("GenBeta.R")
source("PFAdecomp.R")

library(Rlab)
library(MCMCpack)
library(mvtnorm)
library(mnormt)
library(quantreg)
library(robustbase)
library(pracma)

##Part I Simulate beta

P=1000
N=100
NumIter=200

pvec=seq(0.6, 0.9, 0.1) # 

e=0.05
smallestf=0.9

stdized='Y'


betatype="Uninormal"
vartype="v1"

#sigmalist=list("equalcorr4","equalcorr8","FanandSong","IndepCauchy","Nonlinear","Threefactor")
sigmalist=list("IndepCauchy","Nonlinear","Threefactor")

for (mm in 1:length(sigmalist))
{  
 errortype=sigmalist[[mm]]
  
    for(w in pvec) 
    {
      
      z_out_mat=matrix(data=NA, nrow=P, ncol=NumIter)
      var_K_out_mat=matrix(data=NA, nrow=P, ncol=NumIter)
      SD_out_mat=matrix(data=NA, nrow=P, ncol=NumIter)
      hat_beta_out_mat=matrix(data=NA, nrow=P, ncol=NumIter)
      beta_out_mat=matrix(data=NA, nrow=P, ncol=NumIter)
      
      for (iter in 1:NumIter)
      {
       beta=GenBeta(P, w, betatype, vartype)
    
       Sigma=GenSigmaMatrix(N,P,errortype)
    
       hat_beta=mvrnorm(1, beta, Sigma)
       
  
       results_PFA=PFAdecomp(hat_beta, Sigma, smallestf,e,stdized)


       z=results_PFA[[1]]
       
       var_K=results_PFA[[2]]
  
       SD=results_PFA[[3]]

  
       z_out_mat[,iter]=z
       var_K_out_mat[,iter]=var_K 
       SD_out_mat[,iter]=SD
       hat_beta_out_mat[,iter]=hat_beta
       beta_out_mat[,iter]=beta
       
      }
      
      path_name=paste("/./sim_dat_file",betatype,errortype,w,sep="/")
      
      print(path_name)
      file_name=paste(path_name,"z_out_mat.csv",sep="/")
      write.table(z_out_mat, file = file_name, sep = ",", row.names = FALSE,col.names = FALSE )
       
      file_name=paste(path_name,"var_K_out_mat.csv",sep="/")
      write.table(var_K_out_mat, file = file_name, sep = ",", row.names = FALSE,col.names = FALSE )
      
      file_name=paste(path_name,"SD_out_mat.csv",sep="/")
      write.table(SD_out_mat, file = file_name, sep = ",", row.names = FALSE,col.names = FALSE )

      file_name=paste(path_name,"hat_beta_out_mat.csv",sep="/")
      write.table(hat_beta_out_mat, file = file_name, sep = ",", row.names = FALSE,col.names = FALSE )
       
      file_name=paste(path_name,"beta_out_mat.csv",sep="/")
      write.table(beta_out_mat, file = file_name, sep = ",", row.names = FALSE,col.names = FALSE )
    }     
}

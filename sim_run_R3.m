clear;
close all;
clc;

addpath('../LindaZhao_replication/NonParametric');
addpath('../figtree-0.9.1/figtree-0.9.1/matlab');
addpath('../lasso');
addpath('../lasso/sub');

betatype= 'Binormal';   %'Uninormal'; 

sigmatype=  ['Threefactor']; %['Nonlinear'];  %['IndepCauchy']; %['FanandSong'];  %['equalcorr4']; 

delta=[0.1:0.1:0.9]; %0.1;% 0.6

for i=1:length(delta)
   
 z1 = csvread(sprintf('./sim_dat_file/%s/%s/%.1f/z_out_mat.csv',betatype,sigmatype,delta(i)));
 var_k1 = csvread(sprintf('./sim_dat_file/%s/%s/%.1f/var_K_out_mat.csv',betatype,sigmatype,delta(i)));
 hat_beta = csvread(sprintf('./sim_dat_file/%s/%s/%.1f/hat_beta_out_mat.csv',betatype,sigmatype,delta(i)));
 beta = csvread(sprintf('./sim_dat_file/%s/%s/%.1f/beta_out_mat.csv',betatype,sigmatype,delta(i)));
 SD = csvread(sprintf('./sim_dat_file/%s/%s/%.1f/SD_out_mat.csv',betatype,sigmatype,delta(i)));
 
 [n_r, n_c]=size(z1);
 P=n_r;
 T=n_c;
 
 K=0;
 
 for j=1:T
   z_in= z1(:,j);
   var_k1_in=var_k1(:,j);
   hat_beta_in=hat_beta(:,j);
   beta_true=beta(:,j);
   SD_in=SD(:,j);
   
   [w_hat_DepEB_out,mu_hat_DepEB]=EstimationProc(hat_beta_in, z_in, var_k1_in,'DepEB');
   [w_hat_EB_out,mu_hat_EB]=EstimationProc(hat_beta_in, z_in, var_k1_in,'EB');
                                          
  
   mu_hat_DepEB=mu_hat_DepEB'.*SD_in;
   
   w_hat_DepEB(i,j)=w_hat_DepEB_out;
   MSE_DepEB(i,j)=sum((mu_hat_DepEB-beta_true).^2)/P;
   w_hat_EB(i,j)=w_hat_EB_out;
   MSE_EB(i,j)=sum((mu_hat_EB'-beta_true).^2)/P;    
   disp(sprintf('w=%1.2f Trial %d DepEB w=%f EB w=%f',delta(i),j,w_hat_DepEB(i,j),w_hat_EB(i,j)));

 end
 
 
    mean_w_hat_DepEB(i)=mean(w_hat_DepEB(i,:));
    std_w_hat_DepEB(i)=std(w_hat_DepEB(i,:));
    
    mean_w_hat_EB(i)=mean(w_hat_EB(i,:));
    std_w_hat_EB(i)=std(w_hat_EB(i,:));
    
    mean_MSE_DepEB(i)=mean(MSE_DepEB(i,:));
    std_MSE_DepEB(i)=std(MSE_DepEB(i,:));
    
    mean_MSE_EB(i)=mean(MSE_EB(i,:));
    std_MSE_EB(i)=std(MSE_EB(i,:));
 
end

save(sprintf('./sim_dat_file/%s/%s/w_hat_DepEB_%s_%s.mat',betatype,sigmatype,betatype,sigmatype),'w_hat_DepEB');
save(sprintf('./sim_dat_file/%s/%s/w_hat_EB_%s_%s.mat',betatype,sigmatype,betatype,sigmatype),'w_hat_EB');


save(sprintf('sim_dat_file/%s/%s/MSE_DepEB_%s_%s.mat',betatype,sigmatype,betatype,sigmatype),'MSE_DepEB');
save(sprintf('sim_dat_file/%s/%s/MSE_EB_%s_%s.mat',betatype,sigmatype,betatype,sigmatype),'MSE_EB');



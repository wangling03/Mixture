clear;
close all;
clc;

addpath('../matlab_functions/NonParametric');
addpath('../matlab_functions/figtree-0.9.1/figtree-0.9.1/matlab');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tb = readtable('./datafile/data_hc_06.csv');
probenms=Tb(2:end,1);
writetable(probenms,'./datafile/probenames.csv','WriteVariableNames',false) 


beta_PFA = csvread('./datafile/beta_PFA.csv');
var_k_PFA = csvread('./datafile/var_K_PFA.csv');
z = csvread('./datafile/z.csv');
SD = csvread('./datafile/Sigma_SD.csv');



[w_hat_DepEB_out,mu_hat_DepEB]=EstimationProc(z, beta_PFA, var_k_PFA,'DepEB');

numb_sel_DepEB=round((1-w_hat_DepEB_out)*length(beta_PFA)); 
                                                        
                                                         
                                                        
mu_hat_DepEB=mu_hat_DepEB';
[rowsel,colsel]=find(abs(mu_hat_DepEB)>1e-3); 


sel_probe_DepEB=probenms(rowsel,:);

mu_hat_DepEB2=[mu_hat_DepEB, mu_hat_DepEB.*SD];
sel_muhat_DepEB=mu_hat_DepEB2(rowsel,:);


%%%%%%%%%%%%%%%%%%%%%%%%
[w_hat_EB_out,mu_hat_EB]=EstimationProc(z, beta_PFA, var_k_PFA,'EB');

numb_sel_EB=round((1-w_hat_EB_out)*length(z)); %6629 probes

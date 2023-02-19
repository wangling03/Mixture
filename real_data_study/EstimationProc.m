function [w,z_estimate_out]=EstimationProc(beta,z1,var_k1,method)

% z are the observations
% method='DepEB';
% method='EB';
    
%--------------------------------------------------------------------------
% Estimate the hyperparameters
% Compute the estimate given the hyperparameters
%--------------------------------------------------------------------------
    K=0;
    
if strcmp(method,'DepEB')==1

    z=z1./sqrt(var_k1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %scale=std(z);
    scale=median(abs(z)); 
    %scale=1;
    z=z/scale;
    z=full(z);

    [w,g,p_tilde]=NonParametric_EstimateHyperparameters(z);
    [z_estimate]=NonParametric_BayesThreshold_L2_Loss(z,w,g,K);
    z_estimate_out=(scale*z_estimate).*sqrt(var_k1');
end

if strcmp(method,'EB')==1
    
    scale=median(abs(beta)); 
    %scale=1;
    beta=beta/scale;
    beta=full(beta);
    
    [w,g,p_tilde]=NonParametric_EstimateHyperparameters(beta);
    [z_estimate]=NonParametric_BayesThreshold_L2_Loss(beta,w,g,K);
    z_estimate_out=scale*z_estimate;
end

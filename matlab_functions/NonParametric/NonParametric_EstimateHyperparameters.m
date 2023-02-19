function [w,g,p_tilde]=NonParametric_EstimateHyperparameters(x)

%--------------------------------------------------------------------------
% Estimates the hyperparameters by maximizing the marginal likelihood
%--------------------------------------------------------------------------
% INPUTS
%--------------------------------------------------------------------------
% x         ... a vector of n scalar noisy observations
%--------------------------------------------------------------------------
% OUTPUT
%--------------------------------------------------------------------------
% w         ... the mixture parameter of the mixture prior
% g         ... the marginal density of the non-zero part
% p_tilde   ... the posterior probability of being zero
%--------------------------------------------------------------------------
% EXAMPLE USAGE
%--------------------------------------------------------------------------
%  [w,g,p_tilde]=NonParametric_EstimateHyperparameters(x);   
%--------------------------------------------------------------------------
% AUTHORS
%--------------------------------------------------------------------------
% Author: Vikas Chandrakant Raykar and Linda Zhao
% E-Mail: vikas.raykar@siemens.com
% Date  : November 26, 2008
%--------------------------------------------------------------------------

[row,col]=size(x);

if col==1
    x=x';
end

%--------------------------------------------------------------------------
% Optimization paramters
%--------------------------------------------------------------------------

w_tol=1e-6;
options = optimset('Display','off');
iter_max=100;

w_prev=0.5;
p_tilde=0.5*ones(size(x));

w=0.6;

%%figure;

iter=1;
while abs(w-w_prev)>w_tol  & iter < iter_max 
    w_prev=w;
  
    [g]=ComputeMarginalDensity(x,x,p_tilde);
    
    y=[1.5*min(x):0.1:1.5*max(x)];
    [gp]=ComputeMarginalDensity(x,y,p_tilde);
    [gpd]=ComputeMarginalDensityDerivative(x,y,p_tilde);
    %%fig=plot(y,gp,'k');hold on;set(fig,'linewidth',2);
    %%plot(y,gpd,'b');hold on;set(fig,'linewidth',2);
    %%fig=gca; set(fig,'fontsize',14);set(fig,'linewidth',2);
    %%box on;

    w = fminbnd(@(w) ComputeNegativeLogMarginalLikeihood(x,w,g),0,1,options);
    
    
    [p_tilde]=Computeptilde(x,w,g);
    
    disp(sprintf('Iteration %2d w=%1.2f',iter,w));
    %%fig=title(sprintf('Iteration %2d w=%1.2f',iter,w));set(fig,'linewidth',2);
    %print('-dpsc',sprintf('Marginal_iter_%d.eps',iter));
    iter=iter+1;
    %%pause;hold off
    
end

return


function [p_tilde]=Computeptilde(x,w,g)   
%-------------------------------------------------------------------------
% Evaluates p_tilde given x,w, and g
%-------------------------------------------------------------------------

temp=ComputeGaussian(x,0,1);

p_tilde=(w*temp)./((w*temp)+((1-w)*g));

return

function [m]=ComputeNegativeLogMarginalLikeihood(x,w,g)    
%-------------------------------------------------------------------------
% Evaluates the Negative Log Marginal Likelihood of x given
% w and g
%-------------------------------------------------------------------------

m=-sum(log(w*ComputeGaussian(x,0,1)+(1-w)*g));

return

 function [f]=ComputeGaussian(z,mean,variance)
 %-------------------------------------------------------------------------
 % Evaluates the univariate Gaussian with a given mean 
 % and variance.
 %-------------------------------------------------------------------------
 
 f=(1/sqrt(2*pi*variance))*exp(-((z-mean).^2)./(2*variance));
 
 return



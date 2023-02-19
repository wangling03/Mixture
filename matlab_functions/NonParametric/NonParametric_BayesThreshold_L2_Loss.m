function [x]=NonParametric_BayesThreshold_L2_Loss(x,w,g,K)

%--------------------------------------------------------------------------
% Implementation of the thresholding rule for 0/1 + L2 loss
%--------------------------------------------------------------------------
% INPUTS
%--------------------------------------------------------------------------
% x         ... a vector of n scalar noisy observations
% w         ... the mixture parameter of the mixture prior
% g         ... the marginal of the non zero part
% K         ... the penalty term in the mixture-loss function
%--------------------------------------------------------------------------
% OUTPUT
%--------------------------------------------------------------------------
%  x      .... the estimate, thresholded and shrinked observations
%--------------------------------------------------------------------------
% EXAMPLE USAGE
%--------------------------------------------------------------------------
%  [w,g,p_tilde]=NonParametric_EstimateHyperparameters(x); K=0;
%  [x]=NonParametric_BayesThreshold_L2_Loss(x,w,g,K);
%--------------------------------------------------------------------------
% AUTHORS
%--------------------------------------------------------------------------
% Author: Vikas Chandrakant Raykar and Linda Zhao
% E-Mail: vikas.raykar@siemens.com
% Date  : November 27, 2008
%--------------------------------------------------------------------------
% SEE ALSO
%--------------------------------------------------------------------------
% NonParametric_EstimateHyperparameters
%--------------------------------------------------------------------------

[row,col]=size(x);

if col==1
    x=x';
end

[row,col]=size(g);

if col==1
    g=g';
end

%--------------------------------------------------------------------------
% Compute p_tilde
%--------------------------------------------------------------------------

temp=ComputeGaussian(x,0,1);

p_tilde=(w*temp)./((w*temp)+((1-w)*g));

%--------------------------------------------------------------------------
% Compute E
%--------------------------------------------------------------------------

[g_dash]=ComputeMarginalDensityDerivative(x,x,p_tilde);

E=((g_dash)./(g+eps))+x;

%--------------------------------------------------------------------------
% Compute threshold
%--------------------------------------------------------------------------

threshold=(K*(p_tilde)./(((1-p_tilde).^2)));

%--------------------------------------------------------------------------
% Compute the estimate
%--------------------------------------------------------------------------

x=(1-p_tilde).*E;

temp=(E.*E-threshold);

x(temp<0)=0;


return

function [f]=ComputeGaussian(z,mean,variance)

%-------------------------------------------------------------------------
% Evaluates the univariate Gaussian with a given mean
% and variance.
%-------------------------------------------------------------------------

f=(1/sqrt(2*pi*variance))*exp(-((z-mean).^2)./(2*variance));



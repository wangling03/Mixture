function [y]=NonParametric_BayesThreshold_L2_Loss_v2(x,w,g,p_tilde,K,y)

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

[g]=ComputeMarginalDensity(x,y,p_tilde);
[g_dash]=ComputeMarginalDensityDerivative(x,y,p_tilde);

temp=ComputeGaussian(y,0,1);

p_tilde=(w*temp)./((w*temp)+((1-w)*g));

%--------------------------------------------------------------------------
% Compute E
%--------------------------------------------------------------------------


E=y+(g_dash./(g+eps));

%--------------------------------------------------------------------------
% Compute threshold
%--------------------------------------------------------------------------

threshold=(K*(p_tilde)./(((1-p_tilde+eps).^2)));

%--------------------------------------------------------------------------
% Compute the estimate
%--------------------------------------------------------------------------

y=(1-p_tilde).*E;

temp=(E.*E-threshold);

y(temp<0)=0;



return

function [f]=ComputeGaussian(z,mean,variance)

%-------------------------------------------------------------------------
% Evaluates the univariate Gaussian with a given mean
% and variance.
%-------------------------------------------------------------------------

f=(1/sqrt(2*pi*variance))*exp(-((z-mean).^2)./(2*variance));



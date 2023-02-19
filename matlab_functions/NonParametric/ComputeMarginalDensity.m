function [g]=ComputeMarginalDensity(x,y,p_tilde)

%--------------------------------------------------------------------------
% This function computes the weighted kernel estimate of the 
% marginal density of x at y 
% Uses figtree to compute univariate kernel density estimate.
%--------------------------------------------------------------------------

[row,col]=size(x);

if col==1
    x=x';
end

[row,col]=size(y);

if col==1
    y=y';
end

if ~exist('p_tilde','var')
    p_tilde=zeros(size(x));
end

N=length(x);
M=length(y);

% scaling and shifting

min_x=min(x);
min_y=min(y);
shift=min(min_x, min_y); 

x_shifted=x-shift;
y_shifted=y-shift;

max_x=max(x_shifted);
max_y=max(y_shifted);
scale=1/max(max_x,max_y); 

x_shifted_scaled=x_shifted*scale;
y_shifted_scaled=y_shifted*scale;

% Bandwidth selection using the normal reference rule

h=std(x)*((4/(3*N))^(1/5));

%h=std(x)/sqrt(log(N));

%h=std(x)/sqrt(log(N));

%h=0.5;

q=(1/(sqrt(2*pi)*N*h))*(1-p_tilde);

% Fast computation of KDE

epsil=1e-6;

[K,p_max,r]=figtreeChooseParametersNonUniform(1,N,x_shifted_scaled,sqrt(2)*h*scale,epsil,N,1);
[K,rx,clusterIndex,clusterCenter,numPoints,clusterRadii]=figtreeKCenterClustering(1,N,x_shifted_scaled,K);
[p_max]=figtreeChooseTruncationNumber(1,sqrt(2)*h*scale,epsil,rx,1);
[kde]=figtreeEvaluateIfgt(1,N,M,1,x_shifted_scaled,sqrt(2)*h*scale,q,y_shifted_scaled,p_max,K,clusterIndex,clusterCenter,clusterRadii,r,epsil);

g=kde';


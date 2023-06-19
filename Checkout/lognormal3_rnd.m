%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute a random sample of (m * n) values, using the the 3-parameters
% LogNormal distribution. The Hosking and Wallis (1997) version of the
% distribution is chosen. 
%
% Input :
%    -xi (location), alpha (scale), k (shape) : parameters of the 
%    distribution
%    -m, n : number of rows (m) and colums (n) of the output matrix
%
% Output
%   -x : matrix of the random sample
%
% If k > 0, the range of x is : -Inf < x < xi + alpha/k
% If k = 0, the range of x is : -Inf < x < Inf 
% If k < 0, the range of x is : xi + alpha/k <= x < Inf 
%
% If k = 0 , the distribution corresponds to the normal distribution
%
% Source : Hosking, J., & Wallis, J. (1997). Regional Frequency Analysis:
% An Approach Based on L-Moments. Cambridge: Cambridge University Press. 
% doi:10.1017/CBO9780511529443
%
% Guillaume Talbot, INRS-ETE 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=lognormal3_rnd(xi,alpha,k,m,n)
%% Generate a random matrix (m * n) with values between 0 and 1
p=rand(m,n);
%% Compute the inverse cdf of the matrix p to generate the random sample
%Inverse CDF function for the normal distribution
normal_inv=@(P) sqrt(2).*erfinv(2.*P-1);
y=normal_inv(p);
if k==0 %Case Normal distribution
    x=alpha.*y+xi;
else
    x=xi+alpha.*(1-exp(-k.*y))./k;
end
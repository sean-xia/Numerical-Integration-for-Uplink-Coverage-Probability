function pc = UplinkPC_N(t,e,a,l)
%  ref:J. G. Andrews, A. K. Gupta, and H. S. Dhillon, “A Primer on Cellular Network Analysis Using Stochastic Geometry
% input 
%   t: sir threshold
%   e: power control factor, varing from 0 to 1
%   a: path loss index, a>2
%   l：user density
% output
%   coverage probability
if nargin<4
    l = 0.25;  % the default user density
end
tau = 10.^(t/10);
D = @(u,x,r) pi*l.*exp(-l*pi.*u)./(1+tau^(-1).*r.^(a*(e-1)).*u.^(-a*e/2).*x.^a);
% Inner integration of (65)
I1 = @(x,r) integral(@(u) D(u,x,r),0,x.^2,"arrayvalued",true,"AbsTol",1e-6,"RelTol",1e-6);
% vectorization
I1_arr = @(x,r) arrayfun(@(x,r) I1(x,r),x,r);

% The second integration
I2 = @(r) integral(@(x) I1_arr(x,r).*x, 0,20,"arrayvalued",true,"AbsTol",1e-6,"RelTol",1e-6);
% vectorization
I2_arr = @(r) arrayfun(@(r) I2(r),r);

% The outer integration
pc = 2*pi*l*integral(@(r) r.*exp(-pi*l.*r.^2).*exp(-2*pi*l*I2_arr(r)),0,10, ...
    "arrayvalued",true,"AbsTol",1e-6,"RelTol",1e-6);
end
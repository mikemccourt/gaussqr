function [y,M] = ex15b_gqr_truesol(x,t,R,M_sugg)
% function [y,M] = ex15b_gqr_truesol(x,t,R,M_sugg)
% This function evaluates the true solution of the Burgers' equation
% For very large values of R, this uses an asymptotic expansion of the
% modified Bessel function
% 
% Inputs : x - column vector of locations to evaluate at
%          t - time to evaluate at
%          R - viscosity parameter
%          M_sugg - <optional> suggested length of series
% Outputs : y - value of true solution
%           M - length of summation used

x = x(:);
N = length(x);

if nargin==3
    M_sugg = 50;
end

M = M_sugg;
M_mat = repmat(1:M,N,1);
x_mat = repmat(x,1,M);
I_mat = repmat(besseli(1:M,R/(2*pi),1),N,1);
e_mat = exp(-M_mat.^2*pi^2*t/R);
numer = 4*pi/R*sum(M_mat.*I_mat.*sin(pi*x_mat.*M_mat).*e_mat,2);
denom = besseli(0,R/(2*pi),1) + 2*sum(I_mat.*cos(pi*x_mat.*M_mat).*e_mat,2);
y = numer./denom;

% bessasy = @(n,x)(1 - (4*n.^2-1)./(8*x) + ((4*n.^2-1).*(4*n.^2-9))./(128*x.^2) - ((4*n.^2-1).*(4*n.^2-9).*(4*n.^2-25))./(1536*x.^3));
% M = M_sugg;
% M_mat = repmat(1:M,N,1);
% x_mat = repmat(x,1,M);
% I_mat = repmat(bessasy(1:M,R/(2*pi)),N,1);
% e_mat = exp(-M_mat.^2*pi^2*t/R);
% numer = 4*pi/R*sum(M_mat.*I_mat.*sin(pi*x_mat.*M_mat).*e_mat,2);
% denom = bessasy(0,R/(2*pi)) + 2*sum(I_mat.*cos(pi*x_mat.*M_mat).*e_mat,2);
% y = numer./denom;
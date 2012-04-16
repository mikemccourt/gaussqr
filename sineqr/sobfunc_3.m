% Sobolev function
function [MM,s] = sobfunc_3(x,z,L,sigma,beta,t)
% x, z = points (column vector)
% L = value in the interval [0,L]
% sigma = shape parameter
% beta = exponent of the operator
% t = tolerance -> 10^(-t)

sinfunc = @(n,L,x) sqrt(2/L)*sin(pi*x*n/L);

lamfunc = @(n,L,sigma,beta) ((pi*n/L).^2+sigma^2).^(-beta);
 
 M = 1;
 MM = ceil(1/pi*sqrt(10^(t/beta)*((pi*M)^2+(sigma*L)^2)-(sigma*L)^2));
 %ceil(MM), floor(MM)
 
 Xmat = sinfunc(1:MM,L,x);
 Zmat = sinfunc(1:MM,L,z);
 Lmat = diag(lamfunc(1:MM,L,sigma,beta)); size(Lmat);
 Smat = (Xmat.*Zmat)*Lmat;
 
 s = sum(Smat,2);


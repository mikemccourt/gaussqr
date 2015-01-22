% HSSVD_IBBSolve
% This is a compact version of the HSSVD_IBBSolve code for the book

function yeval = HSSVD_IBBSolve(ep,beta,x,y,xeval)
% Evaluate the HS-SVD basis IBB interpolant
% function yeval = HSSVD_IBBSolve_Full(ep,beta,x,y,xeval)
%   Inputs:  ep    - shape parameter
%            beta  - smoothness parameter
%            x     - data locations
%            y     - data values
%            xeval - locations at which to evaluate the interpolant
%   Outputs: yeval - interpolant values at xeval

phifunc = @(n,x) sqrt(2)*sin(pi*x*n);
lamfunc = @(b,e,n) ((pi*n).^2+e^2).^(-b);

N = size(x,1);

M = ceil(1/pi*sqrt(eps^(-1/beta)*(N^2*pi^2+ep^2)-ep^2));
narr = 1:M;

Phi1 = phifunc(narr(1:N),x);
Phi2 = phifunc(narr(N+1:end),x);
lamvec1 = lamfunc(beta,ep,narr(1:N));
lamvec2 = lamfunc(beta,ep,narr(N+1:end));

Rbar = bsxfun(@rdivide,lamvec2',lamvec1).*(Phi2'/Phi1');

Psi = Phi1 + Phi2*Rbar;
c = Psi\y;

Phieval1 = phifunc(narr(1:N),xeval);
Phieval2 = phifunc(narr(N+1:end),xeval);
yeval = Phieval1*c + Phieval2*(Rbar*c);
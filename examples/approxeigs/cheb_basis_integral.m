function Int_Mat = cheb_basis_integral(X,J)
% function Int_Mat = cheb_basis_integral(X,J)
% This function should not be called by the user
% It evaluates the entries in the matrix associated with the Chebyshev
% polynomial basis
%
% NOTE : This does not work yet.  I need to figure out why.

% Evaluate the Chebyshev polynomials
% Don't know how accurate this is for all x in [-1,1]
% Will worry about it at some point
cp = @(n,x) cos(n.*acos(x)).*(n>=0);

% Rescale part of the problem to the orthogonal domain
U = 2*X-1;

Tjm2 = cp(J-2,U);
Tjm1 = cp(J-1,U);
Tjp1 = cp(J+1,U);
Tjp2 = cp(J+2,U);

% Account for possible division by zero
rjm2 = zeros(size(J));
rjm2(J-2>0) = 1./(J(J-2>0)-2);
rjm1 = zeros(size(J));
rjm1(J-1>0) = 1./(J(J-1>0)-1);

rjp1 = 1./(J+1);
rjp2 = 1./(J+2);

m1j   = (-1).^J;
m1jp1 = m1j+1
mm1j  = 1-m1j

Int_Mat = 1/8*( Tjp1.*rjp1 - Tjm1.*rjm1 - 2*m1j.*rjm1.*rjp1 + Tjp2.*rjp2 - Tjm2.*rjm2 - 4*m1j.*rjm2.*rjp2 ) + ...
          1/4*X.*( mm1j.*rjm2.*rjp2 - Tjp1.*rjp1 + Tjm1.*rjm1 + m1jp1.*rjm1.*rjp1 );
      
% We will need to possibly multiply by (-1) for the coefficients to point
% the interpolant in the right direction (positive slope at -1)
% To do this we will exploit Wikipedia
%   Tn'(-1) = (-1)^(n+1)n^2
function Int_Mat = cheb_basis_integral(X,J)
% function Int_Mat = cheb_basis_integral(X,J)
% This function should not be called by the user
% It evaluates the entries in the matrix associated with the Chebyshev
% polynomial basis
%
% We need the J matrix to be ordered [1,2,3,4,...], although at some point
% we could probably remove that restriction

% Evaluate the Chebyshev polynomials
% Don't know how accurate this is for all x in [-1,1]
% Will worry about it at some point
cp = @(n,x) cos(n.*acos(x));

% Rescale part of the problem to the orthogonal domain
U = 2*X-1;

Tjm2 = cp(J-2,U);
Tjm1 = cp(J-1,U);
Tjp1 = cp(J+1,U);
Tjp2 = cp(J+2,U);

% Account for possible division by zero
rjm2 = zeros(size(J));
rjm2(J-2~=0) = 1./(J(J-2~=0)-2);
rjm1 = zeros(size(J));
rjm1(J-1~=0) = 1./(J(J-1~=0)-1);

rjp1 = 1./(J+1);
rjp2 = 1./(J+2);
m1j   = (-1).^J;
mm1j  = 1-m1j;

Int_Mat = 1/16*( 2*Tjp1.*rjp1 - 2*Tjm1.*rjm1 - 4*m1j.*rjm1.*rjp1 + Tjp2.*rjp2 - Tjm2.*rjm2 + 4*m1j.*rjm2.*rjp2 ) + ...
         1/4*X.*(mm1j.*rjm2.*rjp2 - Tjp1.*rjp1 + Tjm1.*rjm1 - mm1j.*rjm1.*rjp1 );
     
% This part of the computation computes the columns associated with the
% first few chebyshev polynomials, which are not acceptable in the
% recursion used above
% In reality, we need to do a search of the J columns, find the ones that
% are associated with the indexes 1,2,3 and do this substitution, but we're
% assuming that the basis functions are ordered [1,2,3,4,...]
IT0 = @(z) z.*(1-z)/2;
IT1 = @(z) (z.^2)/2 -(z.^3)/3-z/6;
IT2 = @(z) -z./6 -(z.^2)/2 + (z.^3)*4/3 -(z.^4)*2/3;
nc = size(Int_Mat,2);
Int_Mat(:,1) = IT0(X(:,1));
if (nc > 1)
    Int_Mat(:,2) = IT1(X(:,2));
    if (nc > 2)
        Int_Mat(:,3) = IT2(X(:,3));
    end
end

     
% We will need to possibly multiply by (-1) for the coefficients to point
% the interpolant in the right direction (positive slope at -1)
% To do this we will exploit Wikipedia
%   Tn'(-1) = (-1)^(n+1)n^2
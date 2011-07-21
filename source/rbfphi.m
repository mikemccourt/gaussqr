% function p = rbfphi(Marr,x,ep,a,b,c,sx2)
% This function should in general be used in the following form
%   function p = rbfphi(Marr,x,ep,a)
% which lets you pass inputs to produce the needed eigenfunctions
%
% Marr - array of exponential values for the eigenfunction
%        [Marr1,Marr2,...,MarrM]
%        Marr1 are the index values for determining the eigenfunctions
%        This should be produced by rbfformMarr
% x - vector of RBF centers
% ep - traditional RBF shape parameter
% a - RBFQR global scale parameter
function p = rbfphi(Marr,x,ep,a,b,c,sx2)

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
logoption = GAUSSQR_PARAMETERS.RBFPHI_WITH_LOGS;

if min(min(Marr))<0
    error(sprintf('Marr=%d is an unacceptable value',min(Marr)))
end
[Mr Mc] = size(Marr);
[xr xc] = size(x);
if Mr~=xc
    error(sprintf('dimension mismatch: Mr=%d, xc=%d',Mr,xc))
end

switch(nargin)
    case {1,2}
        error('Insufficient parameters passed')
    case 3
        a = .5;
        b = ep^2;
        c = sqrt(a^2+2*a*b);
    case 4
        b = ep^2;
        c = sqrt(a^2+2*a*b);
    case 5
        c = sqrt(a^2+2*a*b);
end

if Mc>1
    p = zeros(xr,Mc);
    sx2 = sum(x.^2,2);
    for k=1:Mc
        p(:,k) = rbfphi(Marr(:,k),x,ep,a,b,c,sx2);
    end
else
    if ~exist('sx2') % in case a single n vector is passed
        sx2 = sum(x.^2,2);
    end
    if(b/a<1e-4) % This is a switch for an asymptotic approximation
      q = -.5*(sum(Marr)*log(2)+sum(gammaln(Marr+1))-.25*log(1+2*b/a))+sx2*(-b+.5*b^2/a);
    else
      q = -.5*(sum(Marr)*log(2)+sum(gammaln(Marr+1))-.25*log(1+2*b/a))+sx2*(a-c);
    end
    [Hx,Hi] = HermiteProd(Marr,sqrt(2*c)*x,logoption);
    if logoption
        p = exp(q+Hx).*(-1).^Hi;
    else
        p = exp(q).*Hx;
    end
end
end
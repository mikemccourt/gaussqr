% function p = rbfphialpha(Marr,x,ep,alpha,d,beta,gamma,delta,sx2)
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
% alpha - RBFQR global scale parameter
% d - derivatives in each dimension
% beta,gamma,delta,sx2 - DO NOT PASS, it is for private use
%
% Note: if you pass alpha<0 or ep<0, this will break
function p = rbfphialpha(Marr,x,ep,alpha,d,beta,delta,sx2)

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
logoption = GAUSSQR_PARAMETERS.RBFPHI_WITH_LOGS;
asympttol = GAUSSQR_PARAMETERS.RBFPHI_EXP_TOL;

[Mr Mc] = size(Marr);
[xr xc] = size(x);
if Mr~=xc
    error(sprintf('dimension mismatch: Mr=%d, xc=%d',Mr,xc))
end
if min(min(Marr))<0
    error(sprintf('Marr=%d is an unacceptable value',min(Marr)))
end

switch(nargin)
    case {1,2,3}
        error('Insufficient parameters passed')
    case {4,5}
        if nargin==4
            d = zeros(1,xc); % Default to no derivatives being used
        elseif d>0 || d<0 % Right now you can't use derivatives
            warning(sprintf('%d is an unacceptable derivative, reset to 0',d))
            d = zeros(1,xc);
        end
        if alpha<=0 || ep<=0 || imag(alpha)~=0 || imag(ep)~=0
            error(sprintf('alpha=%g or ep=%g unacceptable; must be real and positive',alpha,ep))
        end
        beta = (1+(2*ep/alpha))^(1/4);
        if beta-1<asympttol % This triggers an asymptotic expansion
            delta = -ep^2+ep^4/alpha^2-2*ep^6/alpha^4;
        else
            delta = 1/2*alpha^2*(1-beta^2);
        end
        p = zeros(xr,Mc);
        sx2 = sum(x.^2,2);
        for k=1:Mc
            p(:,k) = rbfphialpha(Marr(:,k),x,ep,a,d,sx2);
        end
    case {6,7}
        error('Too many parameters passed, see comments')
    case 8
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
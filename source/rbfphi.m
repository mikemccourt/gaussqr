% function p = rbfphi(Marr,x,ep,alpha,deriv,beta,delta,sx2)
% This function should in general be used in the following form
%   function p = rbfphi(Marr,x,ep,alpha)
% which lets you pass inputs to produce the needed eigenfunctions
%
% Marr - array of exponential values for the eigenfunction
%        [Marr1,Marr2,...,MarrM]
%        Marr1 are the index values for determining the eigenfunctions
%        This should be produced by rbfformMarr
% x - vector of RBF centers
% ep - traditional RBF shape parameter
% alpha - RBFQR global scale parameter
% deriv - row vector of derivatives in each dimension
% beta,delta,sx2 - DO NOT PASS, it is for private use
%
% Note: if you pass alpha<0 or ep<0 or m<1, this will error out
%
% Derivatives are calculated using the relationship:
%    d/dx(p_{n}(x)) = 2*delta*x*p_{n}(x)+sqrt(2n-2)*beta*alpha*p_{n-1}(x)
%
% Note: minimum M value is 1, in accordance with new structure
%       anything less than 1 will return 0
function p = rbfphi(Marr,x,ep,alpha,deriv,beta,delta,sx2)

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
logoption = GAUSSQR_PARAMETERS.RBFPHI_WITH_LOGS;
asympttol = GAUSSQR_PARAMETERS.RBFPHI_EXP_TOL;

% Here we define: n as the number of data points
%                 s as the dimension of the data
[Mr Mc] = size(Marr);
[n xc] = size(x);
if Mr~=xc
    error(sprintf('dimension mismatch: Mr=%d, xc=%d',Mr,xc))
else
    s=xc;
end

switch nargin
    case {1,2,3}
        error('Insufficient parameters passed')
    case {4,5}
        if nargin==4
            deriv = zeros(1,s); % Default to no derivatives being used
        else
            [dr dc] = size(deriv);
            if dc~=s
                error(sprintf('dimension mismatch: s=%d, dc=%d',s,dc))
            elseif min(deriv)<0
                error(sprintf('negative derivative requested, %d',min(deriv)))
            elseif dr~=1
                error('d must be a row vector of derivatives')
            end
            if sum(deriv<0 + deriv>2)>0 % Right now you only have 2 derivatives
                warning(sprintf('%d is an unacceptable derivative, reset to 0',deriv))
                deriv = zeros(1,s);
            end
        end
        if alpha<=0 || ep<=0 || imag(alpha)~=0 || imag(ep)~=0
            error(sprintf('alpha=%g or ep=%g unacceptable; must be real and positive',alpha,ep))
        end
        beta = (1+(2*ep/alpha)^2)^(1/4);
        if beta-1<asympttol % This triggers an asymptotic expansion
            delta = -ep^2+ep^4/alpha^2-2*ep^6/alpha^4;
        else
            delta = 1/2*alpha^2*(1-beta^2);
        end
        p = zeros(n,Mc);
        sx2 = sum(x.^2,2);
        for k=1:Mc
            if min(Marr(:,k))<1 % Since H_{-1}=0 by definition
                p(:,k) = 0;
            else
                p(:,k) = rbfphi(Marr(:,k),x,ep,alpha,deriv,beta,delta,sx2);
            end
        end
    case {6,7}
        error('Too many parameters passed, see comments')
    case 8
        if sum(deriv)==0 % No derivatives, easier to handle
            q = -.5*(sum(Marr-1)*log(2)+sum(gammaln(Marr))-s*log(beta))+delta*sx2;
            [Hx,Hi] = HermiteProd(Marr-1,beta*alpha*x,logoption);
            if logoption
                p = exp(q+Hx).*(-1).^Hi;
            else
                p = exp(q).*Hx;
            end
        else
            p = ones(n,1);
            for k=1:s
                m = Marr(k);
                d = deriv(k);
                xk = x(:,k);
                pm = rbfphi(m,xk,ep,alpha,0,beta,delta,sx2);
                switch d % Check which derivative is being asked for
                    case 0
                        p = p.*pm;
                    case 1
                        if m==1
                            pm1 = zeros(size(p));
                        else
                            pm1 = rbfphi(m-1,xk,ep,alpha,0,beta,delta,sx2);
                        end
                        p = p.*(2*delta*x.*pm+sqrt(2*m-2)*beta*alpha*pm1);
                    case 2
                        if m==1
                            pm1 = zeros(size(p));
                            pm2 = zeros(size(p));
                        elseif m==2
                            pm1 = rbfphi(m-1,xk,ep,alpha,0,beta,delta,sx2);
                            pm2 = zeros(size(p));
                        else
                            pm1 = rbfphi(m-1,xk,ep,alpha,0,beta,delta,sx2);
                            pm2 = rbfphi(m-2,xk,ep,alpha,0,beta,delta,sx2);
                        end
                        p = p.*(2*delta*(1+2*delta*xk.^2).*pm + ...
                               4*delta*beta*alpha*sqrt(2*m-2)*xk.*pm1 + ...
                               2*(beta*alpha)^2*sqrt((m-1)*(m-2))*pm2);
                    otherwise
                        error(sprint('Unacceptable derivative %d in rbfphialpha',d))
                end
            end
        end
end

end

% For Developers only
% Warning: This code uses delta<0, as opposed to the analysis which uses
%          delta2>0.  I will change this at some point
% Note: serious improvments may be possible using a recurrence relation -
%       it will generally be the case that we are interested in computing
%       all the eigenfucntions phi_m for 1<m<M, in which case we could
%       exploit the relationship:
%   p_{n+1}(x)=sqrt(2/n)*beta*alpha*x*p_{n}(x)-sqrt(1-1/n)*p_{n-1}(x)
%       Doing so would require us to catch the fact that the user has asked
%       for a sequence of values, and it would also potentially undermine
%       the log thing, since we'd end up computing alot of the values via
%       the recurrence.  The ideal setting would be to compute some values
%       directly and others via the recurrence to balance the savings in
%       cost with the gains from vectorization.
%
%       Of course, along that same thought this whole code could actually
%       probably be improved by computing more than one column of this
%       matrix at a time.  Given that, the dominant cost is still the QR
%       factorization, but it might not always be that way.
%
%       Also, I need to fix this whole issue with the indices being off by
%       one now that we've switched to the alpha evaluation.  I could fix
%       that in the rbfFormMarr function, but I need to think about this a
%       little more.
%
%       Derivatives could potentially be handled combinatorially to allow
%       for more efficient computation, but that would be really hard.  The
%       advantage would be that each term in the summations that represent
%       the derivatives could be handled with logs.
%       Also for derivatives, when we switch to the asymptotic Hermite's
%       there may be a safer way to compute the derivative.  Since the
%       Hermite's have an asymptotic form, that could be differentiated
%       directly rather than using the recurrence relation.
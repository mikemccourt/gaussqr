% function p = rbfphi(Marr,x,ep,alpha,deriv)
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
% alpha - GaussQR global scale parameter
% deriv - row vector of derivatives in each dimension
%         only derivatives of up to 4th order are available
%
% Note: if you pass alpha<0 or ep<0 or m<1, this will error out
%
% Derivatives are calculated using the relationship:
%    d/dx(p_{n}(x)) = -2*delta2*x*p_{n}(x)+sqrt(2n-2)*beta*alpha*p_{n-1}(x)
%
% Note: minimum M value is 1, in accordance with new structure
%       anything less than 1 will return 0
function p = rbfphi(Marr,x,ep,alpha,deriv)

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
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

if nargin<4
    error('Insufficient parameters passed')
elseif nargin==4
    deriv = zeros(1,s); % Default to no derivatives being used
else
    [dr dc] = size(deriv);
    if dc~=s
        error('dimension mismatch: s=%d, dc=%d',s,dc)
    elseif min(deriv)<0
        error('negative derivative requested, %d',min(deriv))
    elseif dr~=1
        error('d must be a row vector of derivatives')
    end
    if sum(deriv<0 + deriv>4)>0 % Right now you only have 4 derivatives
        warning('%d is an unacceptable derivative, reset to 0',deriv)
        deriv = zeros(1,s);
    end
end

if alpha<=0 || ep<=0 || imag(alpha)~=0 || imag(ep)~=0
    error('alpha=%g or ep=%g unacceptable; must be real and positive',alpha,ep)
end

beta = (1+(2*ep/alpha)^2)^(1/4);
if beta-1<asympttol % This triggers an asymptotic expansion
    delta2 = ep^2-ep^4/alpha^2+2*ep^6/alpha^4;
else
    delta2 = 1/2*alpha^2*(beta^2-1);
end

p = zeros(n,Mc);
sx2 = sum(x.^2,2);
for k=1:Mc
    if min(Marr(:,k))>=1 % Since H_{-1}=0 by definition
        p(:,k) = rbfphi_EVAL(Marr(:,k),x,ep,alpha,deriv,beta,delta2,sx2);
    end
end
end






% This below is the private function which actually does the evaluation of
% the eigenfunctions.  You should not ever directly call this.
%
% This function only accepts one multiindex at a time, although I'm sure
% there's a better way to do it.
%
% Derivatives:
%  1 : -2*delta2*xk.*pm+sqrt(2*m-2)*beta*alpha*pm1
%  2 : 2*delta2*(2*delta2*xk.^2-1).*pm + ...
%      -4*sqrt(2)*delta2*beta*alpha*sqrt(m-1)*xk.*pm1 + ...
%      2*(beta*alpha)^2*sqrt((m-1)*(m-2))*pm2;
%  3 : 4*delta2^2*xk.*(3-2*delta2*xk.^2).*pm + ...
%      6*sqrt(2)*delta2*beta*alpha*sqrt(m-1)*(2*delta2*xk.^2-1).*pm1 + ...
%      -12*delta2*(beta*alpha)^2*sqrt((m-1)*(m-2))*xk.*pm2 + ...
%      2*sqrt(2)*(beta*alpha)^3*sqrt((m-1)*(m-2)*(m-3))*pm3
%  4 : 4*delta2^2*(4*delta2^2*xk.^4-12*delta2*xk.^2+3).*pm + ...
%      -16*sqrt(2)*delta2^2*beta*alpha*sqrt(m-1)*xk.*(2*delta2*xk.^2-3).*pm1 + ...
%      24*delta2*(beta*alpha)^2*sqrt((m-1)*(m-2))*(2*delta2*xk.^2-1).*pm2 + ...
%      -16*sqrt(2)*delta2*(beta*alpha)^3*sqrt((m-1)*(m-2)*(m-3))*xk.*pm3 + ...
%      4*(beta*alpha)^4*sqrt((m-1)*(m-2)*(m-3)*(m-4))*pm4
%
% Recurrence Relations:
%  2 : phi_{n-2}(x) = 1/sqrt(m-2) *
%                     (sqrt(2)*beta*alpha*xk .*phi_{n-1}(x) - 
%                      sqrt(m-1)*phi_n(x))
%  3 : phi_{n-3}(x) = 1/(sqrt(m-2)*sqrt(m-3)) *
%                     ((2*(beta*alpha)^2*xk.^2-(n-1)+1) .*phi_{n-1}(x) -
%                      sqrt(2)*sqrt(m-1)*beta*alpha*xk*phi_n(x))
%  4 : phi_{n-4}(x) = 1/(sqrt(m-2)*sqrt(m-3)*(m-4)) *
%                     (sqrt(2)*beta*alpha*xk.*(2*(beta*alpha)^2*xk.^2-2*(m-1)+1) .*phi_{n-1}(x) - 
%                      sqrt(m-1)*(2*(beta*alpha)^2*x.^2-(n-1)+2).*phi_n(x)
%
% The actual derivatives below are computed with a recursion to allow for
% fewer computations of phi
function p = rbfphi_EVAL(Marr,x,ep,a,deriv,b,d2,sx2)
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
logoption = GAUSSQR_PARAMETERS.RBFPHI_WITH_LOGS;

[n s] = size(x);
ba = b*a;
s2 = sqrt(2);

if sum(deriv)==0 % No derivatives, easier to handle
    q = -.5*(sum(Marr-1)*log(2)+sum(gammaln(Marr))-s*log(b))-d2*sx2;
    [Hx,Hi] = HermiteProd(Marr-1,ba*x,logoption);
    if logoption
        p = exp(q+Hx).*(-1).^Hi;
    else
        p = exp(q).*Hx;
    end
else
    if logoption
        Ni = zeros(n,1);
        p = zeros(n,1);
    else
        p = ones(n,1);
    end
    for k=1:s
        m = Marr(k);
        sm1 = sqrt(m-1);
        xk = x(:,k);
        pm = rbfphi_EVAL(m,xk,ep,a,0,b,d2,sx2);
        if m>1
            pm1 = rbfphi_EVAL(m-1,xk,ep,a,0,b,d2,sx2);
        end
        switch deriv(k) % Check which derivative is being asked for
            case 0
                pv = pm;
            case 1
                if m==1
                    pv = -2*d2*xk.*pm;
                else
                    pv = -2*d2*xk.*pm+s2*sm1*ba*pm1;
                end
            case 2
                switch m
                    case 1
                        pv = 2*d2*(2*d2*xk.^2-1).*pm;
                    case 2
                        pv = 2*d2*(2*d2*xk.^2-1).*pm + ...
                             -4*s2*d2*ba*sm1*xk.*pm1;
                    otherwise
                        pv = (2*d2*(2*d2*xk.^2-1)-2*ba^2*(m-1)).*pm + ...
                            (ba^2-2*d2)*2*s2*ba*sm1*xk.*pm1;
                end
            case 3
                switch m
                    case 1
                        pv = 4*d2^2*xk.*(3-2*d2*xk.^2).*pm;
                    case 2
                        pv = 4*d2^2*xk.*(3-2*d2*xk.^2).*pm + ...
                             6*s2*d2*ba*sm1*(2*d2*xk.^2-1).*pm1;
                    case 3
                        pv = -4*d2*xk.*(2*d2^2*x.^2-3*d2-3*ba^2*(m-1)).*pm + ...
                             -6*s2*sm1*d2*ba*(2*xk.^2*ba^2+1-2*d2*xk.^2).*pm1;
                    otherwise
                        pv = (-8*d2^3*xk.^3+(12*d2*ba^2*(m-1)+12*d2^2-4*ba^4*(m-1))*xk).*pm + ...
                             (4*s2*ba*sm1*(3*d2^2+ba^4-3*ba^2*d2)*xk.^2-2*s2*ba*sm1*(ba^2*(m-2)+3*d2)).*pm1;
                end
            case 4
                switch m
                    case 1
                        pv = 4*d2^2*(4*d2^2*xk.^4-12*d2*xk.^2+3).*pm;
                    case 2
                        pv = 4*d2^2*(4*d2^2*xk.^4-12*d2*xk.^2+3).*pm + ...
                             -16*s2*d2^2*ba*sm1*xk.*(2*d2*xk.^2-3).*pm1;
                    case 3
                        pv = 4*d2*(4*d2^3*xk.^4-12*d2*(d2+ba^2*(m-1))*xk.^2+3*(d2+2*ba^2*(m-1))).*pm + ...
                             8*s2*sm1*d2*ba*xk.*(2*d2*(3*ba^2-2*d2)*xk.^2+3*(2*d2-ba^2)).*pm1;
                    case 4
                        pv = 4*d2*(4*d2^3*xk.^4-4*(3*d2^2+3*(m-1)*ba^2*d2-2*ba^4*(m-1))*xk.^2+3*(d2+2*ba^2*(m-1))).*pm + ...
                             -8*s2*sm1*d2*ba*(2*(2*ba^4+2*d2^2-3*d2*ba^2)*xk.^3+(-6*d2+5*ba^2-2*ba^2*(m-1))*xk).*pm1;
                    otherwise
                        pv = 4*(4*d2^4*xk.^4+2*(4*d2*ba^4*(m-1)-6*d2^2*ba^2*(m-1)-ba^6*(m-1)-6*d2^3)*xk.^2+ba^4*(m-3)*(m-1)+3*d2^2+6*d2*ba^2*(m-1)).*pm + ...
                             4*s2*sm1*ba*xk.*(2*(ba^6-4*d2*ba^4+6*d2^2*ba^2-4*d2^3)*xk.^2+4*d2*ba^2*(m-1)+12*d2^2-2*ba^4*(m-1)+3*ba^4-10*d2*ba^2).*pm1;
                end
            otherwise
                error('Unacceptable derivative %d in rbfphi_EVAL, somehow',d)
        end
        if logoption
            Ni = Ni + (pv<0);
            p = p + log(abs(pv));
        else
            p = p.*pv;
        end
    end
    if logoption
        p = exp(p).*(-1).^Ni;
    end
end
end

% For Developers only
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
%           NOTE: I think I fixed this, but I need to confirm
%
%       Derivatives could potentially be handled combinatorially to allow
%       for more efficient computation, but that would be really hard.  The
%       advantage would be that each term in the summations that represent
%       the derivatives could be handled with logs.
%       Also for derivatives, when we switch to the asymptotic Hermite's
%       there may be a safer way to compute the derivative.  Since the
%       Hermite's have an asymptotic form, that could be differentiated
%       directly rather than using the recurrence relation.
%
%       More work should be done to determine the error ramifications of
%       computing derivatives directly (d+1 function evals for d
%       derivatives) versus using recurrence relation
%
%       If logoption is activated, there should be only one step away from
%       the log setting - right at the end.  Everything else should be
%       handled in the log domain.  Right now, 0th derivative terms are
%       computed in the log domain, exponentiated, and then logged again to
%       bring them together in the derivative formulas.  I should make a
%       way to keep the pieces separate and bring them together later.
%
%       Maybe the way to go is to create a function that computes the
%       coefficients for any order derivative in the recurrence form and
%       call that each time it is needed.  Then the eigenfunction and the
%       coefficient computation are separated.

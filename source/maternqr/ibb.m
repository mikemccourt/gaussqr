function [K,M] = ibb(x,z,ep,beta,deriv,Mfix)
% function [K,M] = ibb(x,z,ep,beta,deriv,Mfix)
% This function evaluates the kernel defined by:
%    K_beta(x,z) = sum_1^M lambda_n^beta phi_n(x) phi_n(z)
% where
%    (D^2+ep^2)^beta phi_n = lambda_n phi_n(x)
% using boundary conditions
%    u(0)=u(1)=0, D^2u(0)=D^2u(1)=0, ..., D^{2beta}u(0)=D^{2beta}u(1)=0
% We are calling this function the iterated Brownian Bridge
%
% It should be passed a vector x and an vector z and it will return a
% matrix of size length(x)-by-length(z).  Eventually I will allow this to
% work with higher dimension x and z, but not yet.
%
% When beta=1 or beta=2, we can evaluate a closed form.
% For other parameter values, we can only evaluate it with the summation
% Derivatives are unavailable for beta=1, and only a first derivative is
% available for beta=2.
% Mfix is irrelevant for beta=1 or beta=2.
%
% NOTE: The accuracy of this summation is chosen at 
%          GAUSSQR_PARAMETERS.SUMMATION_TOLERANCE
%       which is defined in rbfsetup.  Choosing a lower tolerance will
%       require more computational time.
%
%%%
%
% Inputs : x - column vector of evaluation points
%          z - column vector of kernel centers
%          ep - peakedness parameter of kernel
%          beta - <optional, default = 1> order of the differential operator
%          deriv - <optional, default = 0> order of derivative
%                  no derivatives greater than 4 can be called
%          Mfix - <optional> demanded length of the summation
% Outputs : K - kernel function evaluation
%           M - <optional> length of summation used
%
%%%
%
% Note that if you pass an x value too long, this probably will crash
% because it will run out of memory.  If that happens, please call it in
% chunks, which will allow it to evaluate successfully.  I may implement a
% try/catch call to allow for this to be executed more gracefully.

if not(exist('deriv','var'))
    deriv = 0;
end

if not(exist('Mfix','var'))
    Mfix = 0;
end

% Make sure all parameters are acceptable
M = IBBKernelErrorCheck(x,z,ep,beta,deriv,Mfix);

if beta==1 && Mfix==0
    mi = bsxfun(@min,x,z');
    ma = bsxfun(@max,x,z')-1;
    K = -sinh(ep*mi).*sinh(ep*ma)/(ep*sinh(1*ep));
elseif beta==2 && Mfix==0
    mi = bsxfun(@min,x,z');
    ma = bsxfun(@max,x,z')-1;
    cem = cosh(ep*mi);
    sem = sinh(ep*mi);
    ceM = cosh(ep*ma);
    seM = sinh(ep*ma);
    if deriv==0
%     K = -1/(2*ep^3*sinh(ep)) * ...
%         (  -ep*min(x,z).*cosh(ep*min(x,z)).*sinh(ep*(max(x,z)-1)) ...
%            -ep*(max(x,z)-1).*sinh(ep*min(x,z)).*cosh(ep*(max(x,z)-1)) ...
%            +(ep*coth(ep)+1)*sinh(ep*min(x,z)).*sinh(ep*(max(x,z)-1))  );
    K = -1/(2*ep^3*sinh(ep)) * ...
        (  -ep*mi.*cem.*seM ...
           -ep*ma.*sem.*ceM ...
           +(ep*coth(ep)+1)*sem.*seM  );
    else
        K = -1/(2*ep^3*sinh(ep)) * (...
            bsxfun(@le,x,z').*(-ep*cem.*seM - ep^2.*mi.*sem.*seM - ep^2*ma.*cem.*ceM + (ep*coth(ep)+1)*ep*cem.*seM ...
                              ) + ...
            bsxfun(@gt,x,z').*(-ep*sem.*ceM - ep^2.*ma.*sem.*seM - ep^2*mi.*cem.*ceM + (ep*coth(ep)+1)*ep*sem.*ceM ...
                              ) ...
                                    );
    end
else
    % Should come up with a way to automate this so that it will check for
    % allocation issues and work around it.
    blocksize = 1000;
    % These are the eigenvalues of the series
    lamfunc = mqr_solveprep();

    rx = size(x,1);
    rz = size(z,1);
    K = zeros(rx,rz);

    % This counts backwards to build up the series kernel from the last
    % term until the n=1 term.
    Mind = M;
    for k=1:floor(M/blocksize)
        n = Mind-blocksize+1:Mind;
        
        Phix = mqr_phi(n,x,1,deriv);
        Phiz = mqr_phi(n,z,1);
        lamvec = lamfunc(n,1,ep,beta);
        K = K + bsxfun(@times,Phix,lamvec)*Phiz';

        Mind = Mind - blocksize;
    end
    if Mind>0
        n = 1:Mind;
        
        Phix = mqr_phi(n,x,1,deriv);
        Phiz = mqr_phi(n,z,1);
        lamvec = lamfunc(n,1,ep,beta);
        K = K + bsxfun(@times,Phix,lamvec)*Phiz';
    end
end



end

function M = IBBKernelErrorCheck(x,z,ep,beta,deriv,Mfix)

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
sumtol = GAUSSQR_PARAMETERS.SUMMATION_TOLERANCE;

[~,cx] = size(x);
[~,cz] = size(z);

if cx~=1 || cz~=1
    error('only column vectors for x and z allowed, size(x,2)=%d, size(z,2)=%d',cx,cz)
elseif ep<0 || length(ep)~=1 || imag(ep)~=0
    error('unacceptable peakedness parameter ep=%g',ep)
elseif ep==0
    error('For ep=0, calls should be made to ppsplinekernel')
end

if min(x)<0 || max(x)>1
    error('This function can only be evaluated on [0,1], min(x)=%g, max(x)=%g',min(x),max(x))
elseif min(z)<0 || max(z)>1
    error('This function must have centers in [0,1], min(z)=%g, max(z)=%g',min(z),max(z))
end

if beta<1 || length(beta)~=1 || imag(beta)~=0 || floor(beta)~=beta
    error('unacceptable order beta=%g',beta)
end
if deriv<0 || deriv>4 || length(deriv)~=1 || imag(deriv)~=0 || floor(deriv)~=deriv
    error('unacceptable derivative deriv=%g',deriv)
end
if Mfix<0 || length(Mfix)~=1 || imag(Mfix)~=0 || floor(Mfix)~=Mfix
    warning('unacceptable length Mfix=%g; pass Mfix=0 to automate decision',Mfix)
end

if deriv>0 && beta==1
    error('First order kernel is not differentiable, deriv=%d, beta=%d',deriv,beta)
elseif deriv>2 && beta==2
    error('Second order kernel can have at most 2 derivatives, deriv=%d, beta=%d',deriv,beta)
end

if Mfix==0
    M = floor(1/pi*sqrt(sumtol^(-1/beta)*(pi^2+ep^2)-ep^2));
else
    M = Mfix;
end

end
function [s,M] = ibb(x,z,ep,beta,deriv,Mfix)
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
sumtol = GAUSSQR_PARAMETERS.SUMMATION_TOLERANCE;

if not(exist('deriv','var'))
    deriv = 0;
end

if not(exist('Mfix','var'))
    Mfix = 0;
end

if Mfix==0
    M = floor(1/pi*sqrt(sumtol^(-1/beta)*(pi^2+ep^2)-ep^2));
else
    M = Mfix;
end

if beta==1 && Mfix==0
    mi = bsxfun(@min,x,z');
    ma = bsxfun(@max,x,z')-1;
    s = -sinh(ep*mi).*sinh(ep*ma)/(ep*sinh(1*ep));
elseif beta==2 && Mfix==0
    mi = bsxfun(@min,x,z');
    ma = bsxfun(@max,x,z')-1;
    cem = cosh(ep*mi);
    sem = sinh(ep*mi);
    ceM = cosh(ep*ma);
    seM = sinh(ep*ma);
    if deriv==0
%     s = -1/(2*ep^3*sinh(ep)) * ...
%         (  -ep*min(x,z).*cosh(ep*min(x,z)).*sinh(ep*(max(x,z)-1)) ...
%            -ep*(max(x,z)-1).*sinh(ep*min(x,z)).*cosh(ep*(max(x,z)-1)) ...
%            +(ep*coth(ep)+1)*sinh(ep*min(x,z)).*sinh(ep*(max(x,z)-1))  );
    s = -1/(2*ep^3*sinh(ep)) * ...
        (  -ep*mi.*cem.*seM ...
           -ep*ma.*sem.*ceM ...
           +(ep*coth(ep)+1)*sem.*seM  );
    else
        s = -1/(2*ep^3*sinh(ep)) * (...
            bsxfun(@le,x,z').*(-ep*cem.*seM - ep^2.*mi.*sem.*seM - ep^2*ma.*cem.*ceM + (ep*coth(ep)+1)*ep*cem.*seM ...
                              ) + ...
            bsxfun(@gt,x,z').*(-ep*sem.*ceM - ep^2.*ma.*sem.*seM - ep^2*mi.*cem.*ceM + (ep*coth(ep)+1)*ep*sem.*ceM ...
                              ) ...
                                    );
    end
else
    % Should come up with a way to automate this so that it will check for
    % allocation issues and work around it.  Maybe it could be optimized to
    % account for cache misses?
    L = 1;
    blocksize = 1000;
    s = zeros(size(x));
    % These are the eigenvalues of the series
    lamfunc = mqr_solveprep();

    rx = size(x,1);

    Mind = M;
    for k=1:floor(M/blocksize)
        n = Mind-blocksize+1:Mind;

        Xmat = mqr_phi(n,x,L,deriv);
        Zmat = mqr_phi(n,z,L);
        Xmat = Xmat.*Zmat;
        Zmat = repmat(lamfunc(n,L,ep,beta),rx,1);
        Xmat = Xmat.*Zmat;

        s = s + sum(Xmat,2);

        Mind = Mind - blocksize;
    end
    if Mind>0
        n = 1:Mind;

        Xmat = mqr_phi(n,x,L,deriv);
        Zmat = mqr_phi(n,z,L);
        Xmat = Xmat.*Zmat;
        Zmat = repmat(lamfunc(n,L,ep,beta),rx,1);
        Xmat = Xmat.*Zmat;

        s = s + sum(Xmat,2);
    end
end



end
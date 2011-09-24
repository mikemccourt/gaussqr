function [alpha,k] = rbfalphasearch(ep,a,b)
% function alpha = rbfalphasearch(ep,a,b)
% ep - shape parameter for the Gaussians
% a - minimum bound(s) of domain
% b - maximum bound(s) of domain
%   In multiple dimensions these should be a=[low_1 ... low_d]
%                                          b=[high_1 ... high_d]
%   These vectors must be the same size, obviously
%
% This function tries to return the smallest alpha for which orthonormality
% is maintained throughout the domain until a set eigenfunction.
%
% You can set the top eigenfunction for which orthonornmality should be
% maintained in rbfsetup with the value
%      GAUSSQR_PARAMETERS.DEFAULT_ORTH_REQUESTED
%
% If this function cannot find an alpha up to what was requested, it will
% try to get as high as possible.  The index up to which orthogonality was
% maintained is returned as k.

% Import global parameters from rbfsetup
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
orthmax = GAUSSQR_PARAMETERS.DEFAULT_ORTH_REQUESTED;
alphamin = GAUSSQR_PARAMETERS.ORTH_MINIMUM_ALPHA;
alphamax = GAUSSQR_PARAMETERS.ORTH_MAXIMUM_ALPHA;
options.TolX = GAUSSQR_PARAMETERS.ORTH_SEARCH_ACCURACY;

k = orthmax;

a = a(:);
b = b(:);

if length(a)~=length(b)
    error('Integration bounds have mismatched dimension: len(a)=%d, len(b)=%d',length(a),length(b))
end

% This is a little dummy iteration to walk through the region which is
% specified by the alphamin and alphamax values.  The first value it finds
% which produces good orthonormality is used as an initial guess for the
% optimization program.  If it gets through the walk without finding a spot
% it will take smaller steps or eventually give up.
goodAreaFound = 0;
intfactor = 5;
while intfactor>1.1 && goodAreaFound==0
    intval = rbforthintegral(k,ep,alphamin,a,b);
    alphaguess = alphamin;
    while intval==0 && alphaguess<alphamax
        alphaguess = alphaguess*intfactor;
        intval = rbforthintegral(k,ep,alphaguess,a,b);
    end
    if alphaguess>alphamax
        intfactor = intfactor/2;
    else
        goodAreaFound = 1;
    end
end

if goodAreaFound == 1
    [alpha,fminval,exitflag] = fminbnd(@(alpha)-rbforthintegral(k,ep,alpha,a,b),alphaguess/intfactor,alphaguess,options);
    if exitflag<1
        warning('Solution not found, consider reducing orthonormality tolerance')
    end
else
    alpha = sqrt(alphamin*alphamax);
    warning('Failed to find an acceptable alpha in [%g,%g]',alphamin,alphamax);
end

% alphavec = logspace(-3,1,50);
% weight = @(a,x) a/sqrt(pi)*exp(-(a*x).^2);
% j = 1;
% erra = zeros(size(alphavec));
% erri = zeros(size(alphavec));
% for al=alphavec
%     erra(j) = quadl(@(x)rbfphialpha(k,x',ep,al)'.^2.*weight(al,x),a(1),b(1))^2;
%     erri(j) = rbforthintegral(k,ep,al,a,b);
%     j = j+1;
% end
% semilogx(alphavec,[erra;-erri])
% [alphavec;erra;erri]

end

function integral = rbforthintegral(k,ep,alpha,a,b)
% Import global parameters from rbfsetup
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
tol = GAUSSQR_PARAMETERS.DEFAULT_ORTH_TOLERANCE;
weight = @(a,x) a/sqrt(pi)*exp(-(a*x).^2);

if alpha<=0 % The optimization could pass whatever
    integral = 0;
else
    intappx = 1;
    if exist('quadgk')
        for j=1:length(a)
            intappx = intappx*quadgk(@(x)rbfphialpha(k,x',ep,alpha)'.^2.*weight(alpha,x),a(j),b(j));
        end
    else
        for j=1:length(a)
            intappx = intappx*quadl(@(x)rbfphialpha(k,x',ep,alpha)'.^2.*weight(alpha,x),a(j),b(j));
        end
    end
    integral = (abs(1-intappx)<tol)*(1/alpha);
end

end

% Developers note: I need to figure out a way to easily test lower function
% indexes if the highest one the user asks for fails
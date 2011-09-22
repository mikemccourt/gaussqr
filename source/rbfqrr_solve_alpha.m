function rbfqrOBJ = rbfqrr_solve(x,y,ep,alpha,M)
% function rbfqrOBJ = rbfqrr_solve(x,y,ep,alpha,M)
% This function accepts required inputs x, y and ep
% x is a Nxd vector of input data points of the form
%     x = [x1';x2';...;xN'], x1 is a column vector
% y is a Nx1 vector of data values at the x locations
% ep is the traditional RBF shape parameter
% alpha is the QR-scale parameter
%
% Optional input is M which is the max length of the
% regression.  If you don't pass M, this will choose an
% M from what is defined in rbfsetup
%
% Note that only positive values of ep and a are used
%
% What is returned from this function is rbfqrOBJ which is
% a large object that encapsulates everything you need to
% do RBF-QR.  It is packaged like this to make it easier for
% you to call other rbfqr functions.
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
Mdefault = GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC;
alphaDefault = GAUSSQR_PARAMETERS.ALPHA_DEFAULT;
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;

if nargin<3
    error('Insufficient inputs')
end
rbfqrOBJ.warnid = '';
rbfqrOBJ.warnmsg = '';

if sum(size(x)~=size(y))
    error('Different sized x (input) and y (output) vectors')
end
if size(x,1)~=size(y,1)
    error('Different numbers of inputs and outputs')
end

if not(exist('alpha'))
    alpha = alphaDefault;
elseif length(alpha)==0
    alpha = alphaDefault;
end

% Checks to make sure that the ep and alpha values are acceptable
if length(ep)>1
    ep = abs(real(ep(1)));
    warning(sprintf('Multiple epsilon values not allowed; using epsilon=%g',ep))
end
if length(alpha)>1
    alpha = abs(real(alpha(1)));
    warning(sprintf('Multiple alpha values not allowed; using alpha=%g',alpha))
end
if abs(real(ep))~=ep
    ep = abs(real(ep));
    warning(sprintf('Only real, positive epsilon allowed; using epsilon=%g',ep))
end
if abs(real(alpha))~=alpha
    alpha = abs(real(alpha));
    warning(sprintf('Only real, positive alpha allowed; using alpha=%g',alpha))
end
if ep==0 || alpha==0
    error(sprintf('Parameters cannot be zero: epsilon=%g, alpha=%g',ep,alpha))
end

[N,d] = size(y);
if Mdefault<=1 % Only a percentage of N is considered
    Mdefault = floor(Mdefault*N);
elseif Mdefault>N
    Mdefault = N;
end

% All this figures out what an acceptable Marr is
% For now I'm doing Marr+1 to solve for an off-by-one issues in the paper
% I'll work on fixing this eventually
if not(exist('M'))
    M = zeros(d,1);
    Mlim = Mdefault;
else
    % Need to check to make sure passed M is reasonable
    if any(M>N)
        warning('M value %d>%d unacceptable, reset to %d',M,N,Mdefault)
        Mlim = Mdefault;
        M = zeros(d,1);
    elseif any(M<0)
        warning('Negative M value passed, setting max Marr length to %d',Mdefault)
        Mlim = Mdefault;
        M = zeros(d,1);
    end
    [Md,Mc] = size(M);
    if Md~=d
        if Md==1 % The user only passed a max length for Marr
            Mlim = M;
            M = zeros(d,1);
        else
            error('M of length %d passed, but data has dimension %d',Md,d)
        end
    elseif Mc~=1
        error('Multiple columns in M detected, but only 1 is allowed')
    else % This is the everything was passed correctly case
        Mlim = N;
    end
end

Marr = rbfformMarr(M,[],Mlim)+1;
phiMat = rbfphialpha(Marr,x,ep,alpha);

lastwarn('')
warning off MATLAB:rankDeficientMatrix
[coef,lsqrrank] = linsolve(phiMat,y);
[warnmsg,msgid] = lastwarn;
if strcmp(msgid,'MATLAB:rankDeficientMatrix')
    rbfqrOBJ.warnid = 'GAUSSQR:lowRankRegression';
    rbfqrOBJ.warnmsg = sprintf('RBFQRr was low rank: rank=%d',lsqrrank);
end
warning on MATLAB:rankDeficientMatrix

rbfqrOBJ.reg   = true;
rbfqrOBJ.ep    = ep;
rbfqrOBJ.alpha = alpha;
rbfqrOBJ.N     = N;
rbfqrOBJ.coef  = coef;
rbfqrOBJ.Marr  = Marr;
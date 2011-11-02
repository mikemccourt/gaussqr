function rbfqrOBJ = rbfqr_solve(x,y,ep,alpha,M)
% function rbfqrOBJ = rbfqr_solve(x,y,ep,alpha,M)
% This function accepts required inputs x, y and ep
% x is a Nxd vector of input data points of the form
%     x = [x1';x2';...;xN'], x1 is a column vector
% y is a Nx1 vector of data values at the x locations
% ep is the traditional RBF shape parameter
%
% There is also a recommended input alpha which is the
% global scale parameter which determines the
% orthogonality enjoyed by the eigenfunctions
% The default is set in rbfsetup, which can be used by passing []
% For the time being alpha is only a single value
% In the future we may allow for different alpha in each dimension
%
% Note that only positive values of ep and alpha are used
%
% Optional input is M>N which is the length of the
% expansion.  If you don't pass M, this will choose an
% M based on the eigenvalues
%
% What is returned from this function is rbfqrOBJ which is
% a large object that encapsulates everything you need to
% do RBF-QR.  It is packaged like this to make it easier for
% you to call other rbfqr functions.

% Import global parameters from rbfsetup
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
Mextramax = GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC;
alphaDefault = GAUSSQR_PARAMETERS.ALPHA_DEFAULT;
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;

if nargin<3
    error('Insufficient inputs')
end
rbfqrOBJ.warnid = '';
rbfqrOBJ.warnmsg = '';

[N,d] = size(x);
if sum(N~=size(y,1))
    error('Different number of input and output points (x and y)')
elseif size(y,2)~=1
    error('You can only pass a 1D output vector y')
end

% If the user didn't pass alpha, or passed [], then we need to pick an
% alpha from rbfalphasearch
computealpha = 0;
if nargin==3
    computealpha = 1;
elseif length(alpha)==0
    computealpha = 1;
end
if computealpha==1
    xminBound = min(x);
    xmaxBound = max(x);
    alpha = rbfalphasearch(ep,xminBound,xmaxBound);
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

nu = (2*ep/alpha)^2;
lam = nu/(2+nu+2*sqrt(1+nu));
if Mextramax<0
    Mextramax = (1-Mextramax/100)*N;
end
MarrN = rbfformMarr(zeros(d,1),[],N);
Mlim = ceil(size(MarrN,2)+log(eps)/log(lam));
% Mlim = ceil(N+log(eps)/log(lam));
if Mextramax==0
    Mextramax = inf; % Allow the array to go as long as it wants
end

% This needs to get better
% Specifically it needs to handle people passing weird stuff
if not(exist('M'))
    M = zeros(d,1);
else
    %M = M(:);
    [Mr Mc] = size(M);
    if Mr~=d
        error('Incorrect M size passed, size(M)=%dx%d d=%d',Mr,Mc,d)
    elseif M<N
        error('rbfqr_solve requires M>N, but M=%g, N=%d',M,N)
    elseif ceil(M)~=M
        warning('Noninteger M passed as %g, reset to %d',M,ceil(M))
        M = ceil(M);
    end
end

Marr = rbfformMarr(M,Mlim,Mextramax)+1;
phiMat = rbfphi(Marr,x,ep,alpha);

[Q,R] = qr(phiMat);
R1 = R(:,1:N);
R2 = R(:,N+1:end);

lastwarn('')
warning off MATLAB:divideByZero
iRdiag = diag(1./diag(R1));
[warnmsg,msgid] = lastwarn;
if strcmp(msgid,'MATLAB:divideByZero')
    rbfqrOBJ.warnid = 'GAUSSQR:zeroQRDiagonal';
    rbfqrOBJ.warnmsg = 'At least one value on the R diagonal was exactly 0';
end
warning on MATLAB:divideByZero

R1s = iRdiag*R1;
opts.UT = true;

lastwarn('')
warning off MATLAB:singularMatrix
if strcmp(rbfqrOBJ.warnid,'GAUSSQR:zeroQRDiagonal')
    Rhat = linsolve(R1s,iRdiag*R2,opts);
else
    Rhat = linsolve(R1s,iRdiag*R2,opts);
    [warnmsg,msgid] = lastwarn;
    if strcmp(msgid,'MATLAB:singularMatrix')
        rbfqrOBJ.warnid = 'GAUSSQR:singularR1invR2';
        rbfqrOBJ.warnmsg = 'Computing inv(R1)R2 ... R1 singular to working precision';
    end
end

D = lam.^(toeplitz(sum(Marr(:,N+1:end),1),sum(Marr(:,N+1:-1:2),1)));
Rbar = D.*Rhat';

[coef,recipcond] = ranksolve(Rhat,Rbar,linsolve(R1s,iRdiag*(Q'*y),opts));
if (recipcond<eps || isnan(recipcond)) && strcmp(rbfqrOBJ.warnid,'')
    rbfqrOBJ.warnid = 'GAUSSQR:illConditionedRanksolve';
    rbfqrOBJ.warnmsg = sprintf('ranksolve encountered an ill-conditioned system, rcond=%g',recipcond);
end
warning on MATLAB:singularMatrix

rbfqrOBJ.reg   = false;
rbfqrOBJ.ep    = ep;
rbfqrOBJ.alpha = alpha;
rbfqrOBJ.N     = N;
rbfqrOBJ.coef  = coef;
rbfqrOBJ.Rbar  = Rbar;
rbfqrOBJ.Marr  = Marr;

if alertuser && ~strcmp(rbfqrOBJ.warnid,'')
    warning(rbfqrOBJ.warnid,rbfqrOBJ.warnmsg)
end

% Developer's note: I should throw something in about handling the
% ill-conditioning.  Users should be alerted when things are bad.
% Specifically it is possible to get some of the diagonal values of R
% exactly equal to zero, which is obviously a problem to which the user
% should be alerted.
% Relevant warning codes include MATLAB:illConditionedMatrix
%                                MATLAB:divideByZero
%                                MATLAB:nearlySingularMatrix
%                                MATLAB:singularMatrix
% I'll need to use warning query to check along the way
%
% Eventually I'll need to find a way to allow users to simultaneously pick
% M and alpha.  Right now I have them guess at M in the rbfsetup file and
% then use that guess to come up with alpha, which then comes up with the
% real M.  This works because there is an acceptable alpha range, but
% there's a better way to do this, I just don't know it yet.
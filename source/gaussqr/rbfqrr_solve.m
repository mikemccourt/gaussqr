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
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;

if nargin<3
    error('Insufficient inputs')
end
rbfqrOBJ.warnid = '';
rbfqrOBJ.warnmsg = '';

[N,d] = size(x);
if N~=size(y,1)
    error('Different numbers of inputs and outputs')
elseif sum(size(y,2)~=1)
    error('Output vector y must have only one column')
end

if nargin==3
    [ep,alpha,Marr] = rbfsolveprep(1,x,ep);
elseif nargin==4
    [ep,alpha,Marr] = rbfsolveprep(1,x,ep,alpha);
else
    [ep,alpha,Marr] = rbfsolveprep(1,x,ep,alpha,M);
end
phiMat = rbfphi(Marr,x,ep,alpha);

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
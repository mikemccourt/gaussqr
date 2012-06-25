function GQR = gqr_rsolve(x,y,ep,alpha,M)
% function GQR = gqr_rolve(x,y,ep,alpha,M)
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
% What is returned from this function is GQR which is
% a large object that encapsulates everything you need to
% do GaussQR.  It is packaged like this to make it easier for
% you to call other rbfqr functions.
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;

% Check input stuff, and call solveprep, to create GQR object
switch nargin
    case {3,4,5}
        [N,d] = size(x);
        if sum(N~=size(y,1))
            error('Different number of input and output points (x and y)')
        elseif size(y,2)~=1
            error('You can only pass a 1D output vector y')
        end
        if nargin==3
            GQR = gqr_solveprep(1,x,ep);
        elseif nargin==4
            GQR = gqr_solveprep(1,x,ep,alpha);
        else
            GQR = gqr_solveprep(1,x,ep,alpha,M);
        end
    otherwise
        error('Unacceptable inputs, nargin=%d',nargin)
end

% Form the linear system
phi = gqr_phi(GQR,x);

% Solve the least squares problem
lastwarn('')
warning off MATLAB:rankDeficientMatrix
[coef,lsqrrank] = linsolve(phi,y);
[warnmsg,msgid] = lastwarn;
if strcmp(msgid,'MATLAB:rankDeficientMatrix')
    GQR.warnid = 'GAUSSQR:lowRankRegression';
    GQR.warnmsg = sprintf('GaussQRr was low rank: rank=%d',lsqrrank);
end
warning on MATLAB:rankDeficientMatrix

% Store the solution
GQR.coef  = coef;

% Tell the user if something went wrong.
if alertuser && ~strcmp(GQR.warnid,'')
    warning(GQR.warnid,GQR.warnmsg)
end
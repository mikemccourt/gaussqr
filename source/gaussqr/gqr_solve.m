function GQR = gqr_solve(x,y,ep,alpha,M)
% function GQR = gqr_solve(x,y,ep,alpha,M)
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
% What is returned from this function is GQR which is
% a large object that encapsulates everything you need to
% do RBF-QR.  It is packaged like this to make it easier for
% you to call other rbfqr functions.

% Import global parameters from rbfsetup
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;
storephi = GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION;

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
            GQR = gqr_solveprep(0,x,ep);
        elseif nargin==4
            GQR = gqr_solveprep(0,x,ep,alpha);
        else
            GQR = gqr_solveprep(0,x,ep,alpha,M);
        end
    otherwise
        error('Unacceptable inputs, nargin=%d',nargin)
end

% Create eigenfunction basis, or recall from previous comptuations
if storephi
	phi = GQR.stored_phi;
else
    phi = gqr_phi(GQR.Marr,x,GQR.ep,GQR.alpha);
end

% Solve in the stable basis
[coef,recipcond] = linsolve(phi*[eye(N);GQR.Rbar],y);
if (recipcond<eps || isnan(recipcond)) && strcmp(GQR.warnid,'')
    GQR.warnid = 'GAUSSQR:illConditionedRanksolve';
    GQR.warnmsg = sprintf('ranksolve encountered an ill-conditioned system, rcond=%g',recipcond);
end

% Store the solution in the GQR object
GQR.coef = coef;

% Tell the user if something went wrong along the way
if alertuser && ~strcmp(GQR.warnid,'')
    warning(GQR.warnid,GQR.warnmsg)
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
%
% I think that the solve may be more stable when the product is computed as
% opposed to it being solved using ranksolve.  Need to check this, and
% reevaluate appropriateness of ranksolve.
%
% Should consider a new approach where the indices are chosen by the rank
% of the R that is created.  Like we should test at each step whether the
% next column is going to add to the rank, or whether it is already
% represented in the space spanned by the earlier columns.  This would
% require some sort of active QR, which checks each column as it comes in.
% I feel like this shouldn't cost more, since all we'd be doing is dumping
% columns that we would otherwise need to be in the solve.
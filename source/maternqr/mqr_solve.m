function MQR = mqr_solve(x,y,L,ep,beta)
% function MQR = mqr_solve(x,y,L,ep,beta)
% This function finds the interpolant to scattered data in 1D using the
% eigenfunction expansion of the Compact Matern kernel
% The problem is defined on the domain [0,L]
%
% Inputs: x - column vector of data points between [0,L]
%         y - column vector of data values evaluated at x
%         L - boundary of the domain, must be positive
%         ep - shape parameter, must be positive
%         beta - smoothness order of spline
% Outputs: MQR - maternQR object with relevant data stored

% Import global parameters from rbfsetup
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;
storephi = GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION;

% Check input stuff, and call solveprep, to create MQR object
switch nargin
    case 5
        [N,d] = size(x);
        if d>1
            error('IBB kernels only work in 1D right now, but size(x)=[%d %d]',size(x))
        elseif sum(N~=size(y,1))
            error('Different number of input and output points (x and y)')
        elseif size(y,2)~=1
            error('You can only pass a 1D output vector y')
        end
        MQR = mqr_solveprep(x,L,ep,beta);
    otherwise
        error('Unacceptable inputs, nargin=%d',nargin)
end

% Create eigenfunction basis, or recall from previous comptuations
if storephi
    Psi = MQR.stored_phi1 + MQR.stored_phi2*MQR.Rbar;
else
    Psi = mqr_phi(MQR,x)*[eye(N);MQR.Rbar];
end

% Solve in the stable basis
[coef,recipcond] = linsolve(Psi,y);
if (recipcond<eps || isnan(recipcond)) && strcmp(MQR.warnid,'')
    MQR.warnid = 'GAUSSQR:illConditionedRanksolve';
    MQR.warnmsg = sprintf('ranksolve encountered an ill-conditioned system, rcond=%g',recipcond);
end

% Store the solution in the MQR object
MQR.coef = coef;

% Tell the user if something went wrong along the way
if alertuser && ~strcmp(MQR.warnid,'')
    warning(MQR.warnid,MQR.warnmsg)
end
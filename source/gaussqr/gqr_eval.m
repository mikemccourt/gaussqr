function [y,GQR] = gqr_eval(GQR,x,deriv)
% function [y,GQR] = gqr_eval(GQR,x,deriv)
% Evaluates a GaussQR approximation at x
%
% Inputs : GQR - GaussQR object created by gqr_solve or gqr_rsolve
%          x - vector of locations at which to evaluate the GQR
%          deriv - <default=0> evaluate the derivative of the GQR
% Outputs : y - value of interpolant at x points
%           GQR - GaussQR object, with stored_phi matrix
%
% Note that if you want to reuse the phi computed here in other evaluations
% that you need to store the GQR object returned here
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;
storephi = GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION;

if ~isfield(GQR,'coef')
    error('no coef field - gqr_eval called before a solve')
end

if nargin==2
    deriv = zeros(1,size(x,2));
end

alreadystored = isfield(GQR,'stored_phi');

if alreadystored && ~storephi
    GQR.warnid = 'GAUSSQR:storageoff';
    GQR.warnmsg = 'Stored x found in GQR, but no storage requested';
    if alertuser
        warning(GQR.warnid,GQR.warnmsg)
    end
end

% Check to see if we need to recompute phi or if we may have already stored
% that matrix earlier
recompute = 0;
if storephi
    if alreadystored % Check if an x was already stored
        if any(size(x)~=size(GQR.stored_x)) || any(size(deriv)~=size(GQR.stored_deriv))
            recompute = 1;
        elseif sum(any(x~=GQR.stored_x)) || any(deriv~=GQR.stored_deriv)
            recompute = 1;
        end
    else % If not, store for the first time
        recompute = 1;
    end
    if recompute
        GQR.stored_x = x;
        GQR.stored_deriv = deriv;
    end
else
    recompute = 1;
end

switch GQR.reg
    case 1
        if recompute
            phiEval = gqr_phi(GQR,x,deriv);
            if storephi
                GQR.stored_phi = phiEval; % Store the phi for later
            end
        else
            phiEval = GQR.stored_phi;
        end
        y = phiEval*GQR.coef;
    case 0
        Rbar = GQR.Rbar;
        N = size(Rbar,2);
        if recompute
            phiEval = gqr_phi(GQR,x,deriv);
            phiEval1 = phiEval(:,1:N);
            phiEval2 = phiEval(:,N+1:end);
            if storephi
                GQR.stored_phi1 = phiEval1; % Store the phi pieces for later
                GQR.stored_phi2 = phiEval2;
            end
        else
            phiEval1 = GQR.stored_phi1;
            phiEval2 = GQR.stored_phi2;
        end
        y = phiEval1*GQR.coef + phiEval2*Rbar*GQR.coef;
    otherwise
        error('reg=%g unacceptable, 0 for interpolation, 1 for regression',reg)
end

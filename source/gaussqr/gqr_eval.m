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

ep    = GQR.ep;
alpha = GQR.alpha;
coef  = GQR.coef;
Marr  = GQR.Marr;
N     = GQR.N;
reg   = GQR.reg;

if nargin==2
    deriv = zeros(1,size(Marr,1));
end

alreadystored = isfield(GQR,'stored_phi');

if alreadystored & ~storephi
    GQR.warnid = 'GAUSSQR:storageoff';
    GQR.warnmsg = 'Stored x found in GQR, but no storage requested';
    if alertuser
        warning(GQR.warnid,GQR.warnmsg)
    end
end

if reg
    if storephi
        if alreadystored % Check if an x was already stored
            if sum(any(x~=GQR.stored_x)) | any(deriv~=GQR.stored_deriv)
                GQR.stored_x = x;
                GQR.stored_deriv = deriv;
                GQR.stored_phi = gqr_phi(Marr,x,ep,alpha,deriv);
            end
        else % If not, store for the first time
            GQR.stored_x = x;
            GQR.stored_deriv = deriv;
            GQR.stored_phi = gqr_phi(Marr,x,ep,alpha,deriv);
        end
        y = GQR.stored_phi*coef;
    else
        phiEval = gqr_phi(Marr,x,ep,alpha,deriv);
        y = phiEval*coef;
    end
else
    Rbar = GQR.Rbar;
    if storephi
        if alreadystored % Check if an x was already stored
            if sum(any(x~=GQR.stored_x)) | any(deriv~=GQR.stored_deriv)
                GQR.stored_x = x;
                GQR.stored_deriv = deriv;
                GQR.stored_phi1 = gqr_phi(Marr(:,1:N),x,ep,alpha,deriv);
                GQR.stored_phi2 = gqr_phi(Marr(:,N+1:end),x,ep,alpha,deriv);
            end
        else % If not, store for the first time
            GQR.stored_x = x;
            GQR.stored_deriv = deriv;
            GQR.stored_phi1 = gqr_phi(Marr(:,1:N),x,ep,alpha,deriv);
            GQR.stored_phi2 = gqr_phi(Marr(:,N+1:end),x,ep,alpha,deriv);
        end
        y = GQR.stored_phi1*coef + GQR.stored_phi2*Rbar*coef;
    else
        phiEval1 = gqr_phi(Marr(:,1:N),x,ep,alpha,deriv);
        phiEval2 = gqr_phi(Marr(:,N+1:end),x,ep,alpha,deriv);
        y = phiEval1*coef + phiEval2*Rbar*coef;
    end
end

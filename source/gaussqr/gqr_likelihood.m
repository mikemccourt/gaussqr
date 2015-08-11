function Lhood = gqr_likelihood(x,y,ep,alpha,lamtrunc)
% This function evaluates the likelihood function using the stable basis.
% Eventually I want this to be part of gqr_solveprep and gqr_solve, but for
% right now it'll just be here.
% Also, for right now this only returns the log.
%
% Realistically, this function is still for only research purposes.
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
storephi = GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION;

if not(exist('lamtrunc'))
    lamtrunc = 0;
end

if ~storephi
    error('Right now you have to store phi for likelihood to work')
end

N = size(x,1);

GQR = gqr_solveprep(0,x,ep,alpha);
phi1 = GQR.stored_phi1;
psi = phi1 + GQR.stored_phi2*GQR.Rbar;


beta = (1+(2*ep/alpha)^2)^.25;
delta2 = alpha^2/2*(beta^2-1);
ead = ep^2 + alpha^2 + delta2;
lamvec = sqrt(alpha^2/ead)*(ep^2/ead).^(0:N-1)';

badlam = find(lamvec<lamtrunc);
laminv = 1./lamvec;
lamvec(badlam) = ones(size(lamvec(badlam)));
laminv(badlam) = zeros(size(laminv(badlam)));
logdetlambda = sum(log(lamvec));


[L,U,P] = lu(phi1);
phi1y = U\(L\(P*y));
logdetphi1 = sum(log(abs(diag(U))));

[L,U,P] = lu(psi);
psiy = U\(L\(P*y));
logdetpsi = sum(log(abs(diag(U))));


GQR_ip = log(phi1y'*(psiy.*laminv));
GQR_det = logdetpsi + logdetphi1 + logdetlambda;
Lhood = -N/2*GQR_ip - 1/2*GQR_det;
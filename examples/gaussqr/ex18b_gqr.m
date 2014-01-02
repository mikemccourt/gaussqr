% ex18b_gqr
% This example tests the validity of the new series expansion approach to
% computing the power function, and compares it to the HS-SVD method and
% the standard basis approach
%
% This does not yet work the way I want it to.  The cancelation still
% occurs and I'm not sure there's any way to avoid it.
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION = 1;

N = 10;
x = pickpoints(-1,1,N);
DM = DistanceMatrix(x,x);

NN = 100;
xx = pickpoints(-1,1,NN);
pnorm = inf;

alpha = 1;

rbf = @(e,r) exp(-(e*r).^2);

Ne = 20;
epvec = logspace(-1,0,Ne);
dvec = zeros(1,Ne);
qvec = zeros(1,Ne);
svec = zeros(1,Ne);

m = 1;
for ep=epvec
    % Compute using the direct method
    K = rbf(ep,DM);
    allk = rbf(ep,DistanceMatrix(xx,x));
    dvec(m) = norm(rbf(ep,0) - sum((allk/K).*allk,2),pnorm);
    
    % Compute using the stably computed cardinal functions
    GQR = gqr_solveprep(0,x,ep,alpha);
    Rbar = GQR.Rbar;
    M = N + size(Rbar,1);
    IRbar = [eye(N);Rbar];
    Phi1 = GQR.stored_phi1;
    Phi2 = GQR.stored_phi2;
    Psi = [Phi1,Phi2]*IRbar;
    allpsi = gqr_phi(GQR,xx)*IRbar;
    qvec(m) = norm(rbf(ep,0) - sum((allpsi/Psi).*allk,2),pnorm);
    
    % Compute using the series expansion
    % Only appropriate for small epsilon
    sLambda1 = sqrt(diag(GQR.eig(1:N)));
    powvals = (allpsi*sLambda1).^2;
    svec(m) = norm(rbf(ep,0) - sum(powvals,2),pnorm);
    
    % This line below lets me test partial sums of the series
%     skvec = ones(NN,N) - cumsum(powvals,2);pause
    
    m = m + 1;
end

% Plot the computed power functions
loglog(epvec,[dvec;qvec;svec])
% ParameterPower
% This example tests different methods for computing the power function of
% a Gaussian interpolation problem.  We use the Gaussian because we have
% access to its HS-SVD

% Define the Gaussian
rbf = @(e,r) exp(-(e*r).^2);

% Create some kernel centers
N = 10;
x = pickpoints(-1,1,N);
Neval = 13;
xeval = pickpoints(-1,1,Neval,'halt');

% Define a region of shape parameters to test
alpha = 1;
epvec = logspace(-1,0,10);
nummat = zeros(Neval,length(epvec));
denvec = zeros(size(epvec));

% Loop through and compute the power function
k = 1;
for ep=epvec
    % Compute the determinant on the interpolation points for epsilon
    GQR = gqr_solveprep(0,x,ep,alpha);
    Phi1 = GQR.stored_phi1;
    Psi = Phi1 + GQR.stored_phi2*GQR.CbarT;
    lamvec = GQR.eig(GQR.Marr(1:N));
    logdetDenom = log(det(Phi1)) + sum(log(lamvec)) + log(det(Psi));
    denvec(k) = logdetDenom;
    
    for m=1:Neval
        GQRx = gqr_solveprep(0,[x;xeval(m)],ep,alpha);
        Phi1x = GQRx.stored_phi1;
        Psix = Phi1x + GQRx.stored_phi2*GQRx.CbarT;
        lamvecx = GQRx.eig(GQRx.Marr(1:N));
        logdetNum = log(det(Phi1x)) + sum(log(lamvecx)) + log(det(Psix));
        nummat(m,k) = real(logdetNum);
    end
    
    k = k + 1;
end
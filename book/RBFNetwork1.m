% RBFNetwork1
% This first example considers Tikhonov Regularization for fixed centers
% and shape parameters
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION = 1;

% Initial example for support-vector machines
if exist('rng','builtin')
    rng(0);
else
    rand('state',0);
    randn('state',0);
end

N = 50;
M = 15;
yf = @(x) (1-4*x+32*x.^2).*exp(-16*x.^2);

% Pick points to evaluate the function at
% Add some error to the data
x = pickpoints(-1,1,N,'rand');
noise = .2;
y = yf(x) + noise*randn(N,1);

% Choose points at which to center the basis functions, as needed
z = pickpoints(-1,1,M,'halton');

% For plotting purposes
xx = pickpoints(-1,1,300);
yy = yf(xx);

% Pick a shape parameter and form the design matrix
ep = .05;
gqr_alpha = 1;
H = rbf(ep,DistanceMatrix(x,z));

% Also, form necessary GQR stuff
% Note that we are feeding this the centers, with which to form the stable
% basis, and also the eigenfunction basis for comparison
GQR = gqr_solveprep(0,z,ep,gqr_alpha);
Phi1 = GQR.stored_phi1;
Lambda1 = diag(GQR.eig(GQR.Marr(1:M)));
Psi = gqr_phi(GQR,x)*[eye(M);GQR.Rbar];
GQR_reg = gqr_solveprep(1,z,ep,gqr_alpha,M);
Phi_reg = gqr_phi(GQR_reg,x);

% Pick a range of Tikhonov Regularization parameters and loop over it
lamvec = logspace(-10,5,40);
dirvec = [];
gqrvec = [];
eigvec = [];
loovec = [];
gcvvec = [];
y_best = 0;
lam_best = 0;
err_best = Inf;
k = 1;
for lam=lamvec
    % Form the variance matrix and solve for the weights
    % This uses the direct method
    A = H'*H + lam*eye(M);
    iAHt = A\(H');
    w = iAHt*y;

    % Evaluate the cost and sum-squared error for that choice
    C = y'*P*y;
    S = y'*P*P*y;

    % Evaluate predictions on the test/plotting points
    % Check the error
    yp = rbf(ep,DistanceMatrix(xx,z))*w;
    dirvec(k) = errcompute(yp,yy);

    % Evaluate the parameterization schemes
    % The projection matrix is needed for this
    P = eye(N) - H*iAHt;
    loovec(k) = y'*P*diag(1./diag(P).^2)*P*y/N;
    gcvvec(k) = N*y'*P^2*y/trace(P)^2;
    
    % Compute instead the coefficients for the HS-SVD method
    % We will first compute with the stable basis
    GQR.coef = (Psi'*Psi + lam*eye(M))\(Psi'*y);
    yp = gqr_eval(GQR,xx);
    gqrvec(k) = errcompute(yp,yy);
    
    % Compute as well with the eigenfunctions
    GQR_reg.coef = (Phi_reg'*Phi_reg + lam*eye(M))\(Phi_reg'*y);
    yp = gqr_eval(GQR_reg,xx);
    eigvec(k) = errcompute(yp,yy);
    
    if gqrvec(k)<err_best
        err_best = gqrvec(k);
        lam_best = lam;
        y_best = yp;
    end
    k = k + 1;
end
loglog(lamvec,[dirvec;gqrvec;eigvec;gcvvec])
title(sprintf('ep=%g,N=%d,M=%d',ep,N,M))
legend('Direct','Stable','Eigs','GCV')

GQR_reg.coef = Phi_reg\y;
yp = gqr_eval(GQR_reg,xx);
lam0_err = errcompute(yp,yy);

% Plot the results
figure
plot(x,y,'or')
hold on
plot(xx,yy,'linewidth',2)
plot(xx,y_best,'--k','linewidth',2)
plot(xx,yp,'-.m','linewidth',2)
hold off
ylim([-1,2])
title(sprintf('ep=%g,lam=%g,N=%d,M=%d',ep,lam_best,N,M))
legend('Data','True',sprintf('lam=%g err=%g',lam_best,err_best),sprintf('lam=0 err=%g',lam0_err))
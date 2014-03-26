% RBFNetwork2
% This example considers a fixed number of basis functions and studies the
% effect of varying ep and lam values
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
M = 10;
yf = @(x) (1-4*x+32*x.^2).*exp(-16*x.^2);

% Pick points to evaluate the function at
% Add some error to the data
x = pickpoints(-1,1,N,'rand');
noise = .2;
y = yf(x) + noise*randn(N,1);

% Choose points at which to center the basis functions
z = pickpoints(-1,1,M);

% For plotting purposes
xx = pickpoints(-1,1,300);
yy = yf(xx);

% Pick a range of shape and regularization parameters
gqr_alpha = 1;
N_ep = 30;
N_lam = 35;
epvec = logspace(-2,1,N_ep);
lamvec = logspace(-10,5,N_lam);

dirmat = zeros(N_ep,N_lam);
gqrmat = zeros(N_ep,N_lam);
eigmat = zeros(N_ep,N_lam);
err_best = Inf;
k = 1;
for ep=epvec
    % Pick a shape parameter and form the design matrix
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
    
    j = 1;
    for lam=lamvec
        % Form the variance matrix and solve for the weights
        % This uses the direct method
        A = H'*H + lam*eye(M);
        w = A\(H'*y);

        % Evaluate predictions on the test/plotting points
        % Check the error
        yp = rbf(ep,DistanceMatrix(xx,z))*w;
        dirmat(k,j) = errcompute(yp,yy);

        % Compute instead the coefficients for the HS-SVD method
        % We will first compute with the stable basis
        GQR.coef = (Psi'*Psi + lam*eye(M))\(Psi'*y);
        yp = gqr_eval(GQR,xx);
        gqrmat(k,j) = errcompute(yp,yy);

        % Compute as well with the eigenfunctions
        GQR_reg.coef = (Phi_reg'*Phi_reg + lam*eye(M))\(Phi_reg'*y);
        yp = gqr_eval(GQR_reg,xx);
        eigmat(k,j) = errcompute(yp,yy);

        if gqrmat(k,j)<err_best
            err_best = gqrmat(k,j);
            ep_best = ep;
            lam_best = lam;
            y_best = yp;
        end
        j = j + 1;
    end
    k = k + 1
end
[E,L] = meshgrid(epvec,lamvec);
subplot(1,3,1)
h = surf(E,L,log10(dirmat'));
title(sprintf('standard,N=%d,M=%d',N,M))
set(h,'edgecolor','none')
set(gca,'xscale','log');set(gca,'yscale','log');
xlim([1e-2,1e1]),ylim([1e-10,1e5]),zlim([-2.5,-1])
xlabel('\epsilon'),ylabel('\lambda'),zlabel('log_{10} error')
subplot(1,3,2)
h = surf(E,L,log10(gqrmat'));
title(sprintf('stable,N=%d,M=%d',N,M))
set(h,'edgecolor','none')
set(gca,'xscale','log');set(gca,'yscale','log');
xlim([1e-2,1e1]),ylim([1e-10,1e5]),zlim([-2.5,-1])
xlabel('\epsilon'),ylabel('\lambda'),zlabel('log_{10} error')
subplot(1,3,3)
h = surf(E,L,log10(eigmat'));
title(sprintf('eigs,N=%d,M=%d',N,M))
set(h,'edgecolor','none')
set(gca,'xscale','log');set(gca,'yscale','log');
xlim([1e-2,1e1]),ylim([1e-10,1e5]),zlim([-2.5,-1])
xlabel('\epsilon'),ylabel('\lambda'),zlabel('log_{10} error')

GQR_reg = gqr_solveprep(1,z,1e-8,gqr_alpha,M);
Phi_reg = gqr_phi(GQR_reg,x);
GQR_reg.coef = Phi_reg\y;
yp = gqr_eval(GQR_reg,xx);
lamep0_err = errcompute(yp,yy);

% Plot the results
figure
plot(x,y,'or')
hold on
plot(xx,yy,'linewidth',2)
plot(xx,y_best,'--k','linewidth',2)
plot(xx,yp,'-.m','linewidth',2)
hold off
ylim([-1,2])
title(sprintf('N=%d,M=%d',N,M))
legend('Data','True',sprintf('lam=%g ep=%g err=%g',lam_best,ep_best,err_best),sprintf('lam=ep=0 err=%g',lamep0_err))
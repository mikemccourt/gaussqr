% RBFNetwork2
% This example considers a fixed number of basis functions and studies the
% effect of varying ep and lam values
% Specifically, this is studying the GCV value

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
x = pickpoints(-1,1,N-2,'rand');
x = [x;-1;1];
noise = .2;
y = yf(x) + noise*randn(N,1);

% Choose points at which to center the basis functions
z = pickpoints(-1,1,M);

% For plotting purposes
xx = pickpoints(-1,1,300);
yy = yf(xx);

% Pick a range of shape and regularization parameters
gqr_alpha = 1;
N_ep = 50;
N_lam = 55;
epvec = logspace(-2,1,N_ep);
lamvec = logspace(-10,5,N_lam);

loomat = zeros(N_ep,N_lam);
gcvmat = zeros(N_ep,N_lam);
eigmat = zeros(N_ep,N_lam);
err_eig = Inf;
err_gcv = Inf;
err_loo = Inf;
k = 1;
for ep=epvec
    % Also, form necessary GQR stuff
    % Note that we are using the eigenfunctions only here
    GQR = gqr_solveprep(1,z,ep,gqr_alpha,M);
    Phi = gqr_phi(GQR,x);
    
    j = 1;
    for lam=lamvec
        % Form the variance matrix and solve for the weights
        % This uses the direct method
        H = Phi;
        A = H'*H + lam*eye(M);
        iAHt = A\(H');

        % Evaluate the cost and sum-squared error for that choice
        P = eye(N) - H*iAHt;
        Py = P*y;
        C = y'*Py;
        S = Py'*Py;

        % Evaluate the parameterization schemes
        % The projection matrix is needed for this
        loomat(k,j) = Py'*diag(1./diag(P).^2)*Py/N;
        gcvmat(k,j) = N*S/trace(P)^2;

        % Evaluate predictions on the test/plotting points
        % Check the error
        GQR.coef = iAHt*y;
        yp = gqr_eval(GQR,xx);
        eigmat(k,j) = errcompute(yp,yy);

        if eigmat(k,j)<err_eig
            err_eig = eigmat(k,j);ep_eig = ep;lam_eig = lam;y_eig = yp;
        end
        if gcvmat(k,j)<err_gcv
            err_gcv = gcvmat(k,j);ep_gcv = ep;lam_gcv = lam;y_gcv = yp;
        end
        if loomat(k,j)<err_loo
            err_loo = loomat(k,j);ep_loo = ep;lam_loo = lam;y_loo = yp;
        end
        j = j + 1;
    end
    k = k + 1
end

% Plot the errors of the various solution strategies
[E,L] = meshgrid(epvec,lamvec);
subplot(1,3,1)
h = surf(E,L,log10(eigmat'));
title(sprintf('Error,N=%d,M=%d',N,M))
set(h,'edgecolor','none')
set(gca,'xscale','log');set(gca,'yscale','log');
xlim([1e-2,1e1]),ylim([1e-10,1e5]),zlim([-2.5,-1])
xlabel('\epsilon'),ylabel('\lambda'),zlabel('log_{10} error')
subplot(1,3,2)
h = surf(E,L,log10(loomat'));
title(sprintf('LOOCV,N=%d,M=%d',N,M))
set(h,'edgecolor','none')
set(gca,'xscale','log');set(gca,'yscale','log');
xlim([1e-2,1e1]),ylim([1e-10,1e5])
xlabel('\epsilon'),ylabel('\lambda'),zlabel('log_{10} loocv')
subplot(1,3,3)
h = surf(E,L,log10(gcvmat'));
title(sprintf('GCV,N=%d,M=%d',N,M))
set(h,'edgecolor','none')
set(gca,'xscale','log');set(gca,'yscale','log');
xlim([1e-2,1e1]),ylim([1e-10,1e5])
xlabel('\epsilon'),ylabel('\lambda'),zlabel('log_{10} gcv')

% Also solve the problem in the flat limit, with no regularization
GQR = gqr_solveprep(1,z,1e-8,gqr_alpha,M);
Phi = gqr_phi(GQR,x);
GQR.coef = Phi\y;
yp = gqr_eval(GQR,xx);
lamep0_err = errcompute(yp,yy);

% Plot the results
figure
plot(x,y,'or')
hold on
plot(xx,yy,'linewidth',2)
plot(xx,y_eig,'--k','linewidth',2)
plot(xx,y_loo,'--c','linewidth',2)
plot(xx,y_gcv,'--g','linewidth',2)
plot(xx,yp,'-.m','linewidth',2)
hold off
ylim([-1,2])
title(sprintf('N=%d,M=%d',N,M))
legend('Data','True',...
    sprintf('best lam=%4g ep=%4g err=%4g',lam_eig,ep_eig,err_eig),...
    sprintf('LOO lam=%4g ep=%4g err=%4g', lam_loo,ep_loo,err_loo),...
    sprintf('GCV lam=%4g ep=%4g err=%4g', lam_gcv,ep_gcv,err_gcv),...
    sprintf('lam=ep=0 err=%g',lamep0_err))
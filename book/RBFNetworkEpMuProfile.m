% RBFNetworkEpMuProfile
% We study the effect of different parameterization techniques over a two
% dimensional region of epsilon and mu values

global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.RANDOM_SEED(0);

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
N_mu = 55;
epvec = logspace(-2,1,N_ep);
muvec = logspace(-10,5,N_mu);

loomat = zeros(N_ep,N_mu);
gcvmat = zeros(N_ep,N_mu);
eigmat = zeros(N_ep,N_mu);
err_eig = Inf;
gcv_min = Inf;
loo_min = Inf;
k = 1;
h_waitbar = waitbar(0,'Initiating');pause(.1)
progress = 0;
for ep=epvec
    % Form the Gaussian eigenfunction object and evaluate/decompose Phi
    GQR = gqr_solveprep(1,z,ep,gqr_alpha,M);
    Phi = gqr_phi(GQR,x);
    [U,S,V] = svd(Phi,0);Sv = diag(S);
    
    j = 1;
    for mu=muvec
        % Evaluate the projection matrix and residual for the CV calculations
        P = eye(N) - U*diag(1./(1+mu./Sv.^2))*U';
        Py = P*y;

        % Evaluate the CV parameterization schemes
        loomat(k,j) = Py'*diag(1./diag(P).^2)*Py/N;
        gcvmat(k,j) = N*(Py'*Py)/trace(P)^2;

        % Evaluate predictions on the test/plotting points
        % Check the error
        GQR.coef = V*((U'*y)./(Sv+mu./Sv));
        yp = gqr_eval(GQR,xx);
        eigmat(k,j) = errcompute(yp,yy);

        % Record any optimal values that occur during the search
        if eigmat(k,j)<err_eig
            err_eig = eigmat(k,j);ep_eig = ep;mu_eig = mu;y_eig = yp;
        end
        if gcvmat(k,j)<gcv_min
            gcv_min = gcvmat(k,j);ep_gcv = ep;mu_gcv = mu;y_gcv = yp;
        end
        if loomat(k,j)<loo_min
            loo_min = loomat(k,j);ep_loo = ep;mu_loo = mu;y_loo = yp;
        end
        j = j + 1;
    end
    k = k + 1;
    progress = floor(100*k/N_ep)/100;
    waitbar(progress,h_waitbar,'Computing')
end
waitbar(1,h_waitbar,'Plotting')

% Plot the errors of the various solution strategies
h_subplots = figure;
[E,L] = meshgrid(epvec,muvec);
% This is the plot of the computed error
subplot(1,3,1)
h = surf(E,L,log10(eigmat'));
% title(sprintf('Error,N=%d,M=%d',N,M))
set(h,'edgecolor','none')
set(gca,'xscale','log');set(gca,'yscale','log');
xlim([1e-2,1e1]),ylim([1e-10,1e5]),zlim([-2.5,-1])
xlabel('\epsilon'),ylabel('\mu'),zlabel('log_{10} error')
view([-1 -1 1]),grid off
% This is the plot of the LOOCV valuess
subplot(1,3,2)
h = surf(E,L,log10(loomat'));
% title(sprintf('LOOCV,N=%d,M=%d',N,M))
set(h,'edgecolor','none')
set(gca,'xscale','log');set(gca,'yscale','log');
xlim([1e-2,1e1]),ylim([1e-10,1e5])
xlabel('\epsilon'),ylabel('\mu'),zlabel('log_{10} loocv')
view([-1 -1 1]),grid off
% This is the plot of the GCV values
subplot(1,3,3)
h = surf(E,L,log10(gcvmat'));
% title(sprintf('GCV,N=%d,M=%d',N,M))
set(h,'edgecolor','none')
set(gca,'xscale','log');set(gca,'yscale','log');
xlim([1e-2,1e1]),ylim([1e-10,1e5])
xlabel('\epsilon'),ylabel('\mu'),zlabel('log_{10} gcv')
view([-1 -1 1]),grid off
% This will apply to all three graphs
% Pink makes it easier to see difference near bottom
% Also changing the background color to not white
colormap('pink'),whitebg(h_subplots,[.7 .8 .8])

% Also solve the problem in the flat limit, with no regularization
GQR = gqr_solveprep(1,z,1e-8,gqr_alpha,M);
Phi = gqr_phi(GQR,x);
GQR.coef = Phi\y;
yp = gqr_eval(GQR,xx);

% Compute the necessary errors (not the mins of GCV, LOOCV)
err_muep0 = errcompute(yp,yy);
err_gcv = errcompute(y_gcv,yy);
err_loo = errcompute(y_loo,yy);

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
    sprintf('best mu=%2.2g ep=%2.2g err=%2.2g',mu_eig,ep_eig,err_eig),...
    sprintf('LOO mu=%2.2g ep=%2.2g err=%2.2g', mu_loo,ep_loo,err_loo),...
    sprintf('GCV mu=%2.2g ep=%2.2g err=%2.2g', mu_gcv,ep_gcv,err_gcv),...
    sprintf('mu=ep=0 err=%2.2g',err_muep0))

close(h_waitbar)
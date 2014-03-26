% RBFNetwork2
% This example considers a fixed number of basis functions and studies the
% effect of varying ep and lam values

% Initial example for support-vector machines
if exist('rng','builtin')
    rng(0);
else
    rand('state',0);
    randn('state',0);
end

N = 50;
yf = @(x) (1-4*x+32*x.^2).*exp(-16*x.^2);
rbf = @(e,r) exp(-(e*r).^2);

% Pick points to evaluate the function at
% Add some error to the data
% Make sure the end points are included because we can't interpolate beyond
% the end points - that would be extrapolation
x = pickpoints(-1,1,N-2,'rand');
x = [x;-1;1];
noise = .2;
y = yf(x) + noise*randn(N,1);

% For plotting purposes
xx = pickpoints(-1,1,300);
yy = yf(xx);

% Pick a range of shape and regularization parameters
gqr_alpha = 1;
N_ep = 30;
N_lam = 35;
epvec = logspace(-2,1,N_ep);
lamvec = logspace(-10,5,N_lam);
Mvec = [10,15,20];

m = 1;
err_best = Inf;
for M=Mvec
    % Choose points at which to center the basis functions
    z = pickpoints(-1,1,M);
    
    dirmat = zeros(N_ep,N_lam);
    gqrmat = zeros(N_ep,N_lam);
    eigmat = zeros(N_ep,N_lam);
    k = 1;
    for ep=epvec
        % Pick a shape parameter and form the design matrix
        H = rbf(ep,DistanceMatrix(x,z));

        % Also, form necessary GQR stuff
        % Note that we are feeding this the centers, with which to form the stable
        % basis, and also the eigenfunction basis for comparison
        GQR = gqr_solveprep(1,z,ep,gqr_alpha,M);
        Phi = gqr_phi(GQR,x);

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

            % Compute as well with the eigenfunctions
            GQR.coef = (Phi'*Phi + lam*eye(M))\(Phi'*y);
            yp = gqr_eval(GQR,xx);
            eigmat(k,j) = errcompute(yp,yy);

            if eigmat(k,j)<err_best
                err_best = eigmat(k,j);
                M_best = M;
                ep_best = ep;
                lam_best = lam;
                y_best = yp;
            end
            j = j + 1;
        end
        k = k + 1;
    end

    % Plot the errors of the various solution strategies
    [E,L] = meshgrid(epvec,lamvec);
    subplot(2,3,m)
    h = surf(E,L,log10(dirmat'));
    title(sprintf('standard,N=%d,M=%d',N,M))
    set(h,'edgecolor','none')
    set(gca,'xscale','log');set(gca,'yscale','log');
    xlim([1e-2,1e1]),ylim([1e-10,1e5]),zlim([-2.5,-1])
    xlabel('\epsilon'),ylabel('\lambda'),zlabel('log_{10} error')
    subplot(2,3,m+length(Mvec))
    h = surf(E,L,log10(eigmat'));
    title(sprintf('eigs,N=%d,M=%d',N,M))
    set(h,'edgecolor','none')
    set(gca,'xscale','log');set(gca,'yscale','log');
    xlim([1e-2,1e1]),ylim([1e-10,1e5]),zlim([-2.5,-1])
    xlabel('\epsilon'),ylabel('\lambda'),zlabel('log_{10} error')
    
    m = m + 1
end

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
plot(xx,y_best,'--k','linewidth',2)
plot(xx,yp,'-.m','linewidth',2)
hold off
ylim([-1,2])
title(sprintf('N=%d,M=%d',N,M))
legend('Data','True',sprintf('M=%d lam=%g ep=%g err=%g',M_best,lam_best,ep_best,err_best),sprintf('M=%d lam=ep=0 err=%g',M,lamep0_err))
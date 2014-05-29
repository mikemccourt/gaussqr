% RBFNetwork1
% This first example considers Tikhonov Regularization for fixed centers
% and shape parameters

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
rbf = @(e,r) exp(-(e*r).^2);

% Pick points to evaluate the function at
% Add some error to the data
x = pickpoints(-1,1,N-2,'rand');
x = [x;-1;1];
noise = .2;
y = yf(x) + noise*randn(N,1);

% Choose points at which to center the basis functions, as needed
z = pickpoints(-1,1,M,'halton');

% For plotting purposes
xx = pickpoints(-1,1,300);
yy = yf(xx);

% Pick a shape parameter and form the design matrix
ep = .1;
gqr_alpha = 1;
H = rbf(ep,DistanceMatrix(x,z));

% Also, form necessary GQR stuff
% Note that we are feeding this the centers, with which to form the stable
% basis, and also the eigenfunction basis for comparison
GQR = gqr_solveprep(0,z,ep,gqr_alpha);
Psi = gqr_phi(GQR,x)*[eye(M);GQR.Rbar];
GQR_reg = gqr_solveprep(1,z,ep,gqr_alpha,M);
Phi_reg = gqr_phi(GQR_reg,x);

% Pick a range of Tikhonov Regularization parameters and loop over it
lamvec = logspace(-15,5,40);
dirvec = [];gqrvec = [];eigvec = [];
loovec = [];loevec = [];
gcvvec = [];gcevec = [];
eig_err_best = Inf;dir_err_best = Inf;gce_err_best = Inf;gcd_err_best = [];
k = 1;
for lam=lamvec
    % Form the variance matrix and solve for the weights
    % This uses the direct method
    A = H'*H + lam*eye(M);
    iAHt = A\(H');
    w = iAHt*y;

    % Evaluate the cost and sum-squared error for that choice
    P = eye(N) - H*iAHt;
    Py = P*y;
    S = Py'*Py;

    % Evaluate predictions on the test/plotting points
    % Check the error
    yp = rbf(ep,DistanceMatrix(xx,z))*w;
    dirvec(k) = errcompute(yp,yy);

    % Evaluate the parameterization schemes
    % The projection matrix is needed for this
    lodvec(k) = Py'*diag(1./diag(P).^2)*Py/N;
    gcdvec(k) = N*S/trace(P)^2;
    
    % Compute instead the coefficients for the HS-SVD method
    % We will first compute with the stable basis
    GQR.coef = (Psi'*Psi + lam*eye(M))\(Psi'*y);
    yp = gqr_eval(GQR,xx);
    gqrvec(k) = errcompute(yp,yy);
    
    % Compute as well with the eigenfunctions
    % First the CV content for the eigenfunction basis
    % Then the predictions at the test points
    A = (Phi_reg'*Phi_reg + lam*eye(M));
    iAHt = A\(Phi_reg');
    w = iAHt*y;
    P = eye(N) - Phi_reg*iAHt;
    Py = P*y;
    gcevec(k) = N*Py'*Py/trace(P)^2;
    loevec(k) = Py'*diag(1./diag(P).^2)*Py/N;
    GQR_reg.coef = w;
    yp = gqr_eval(GQR_reg,xx);
    eigvec(k) = errcompute(yp,yy);

    if eigvec(k)<eig_err_best
        eig_err_best = eigvec(k);eig_lam_best = lam;eig_y_best = yp;
    end
    if dirvec(k)<dir_err_best
        dir_err_best = dirvec(k);dir_lam_best = lam;dir_y_best = yp;
    end
    if gcevec(k)<gce_err_best
        gce_err_best = gcevec(k);gce_lam_best = lam;gce_y_best = yp;
    end
    if gcdvec(k)<gce_err_best
        gcd_err_best = gcdvec(k);gcd_lam_best = lam;gcd_y_best = yp;
    end
    k = k + 1;
end

[tmp,id] = min(dirvec);
[tmp,ig] = min(gcdvec);
[tmp,ie] = min(eigvec);
[tmp,ic] = min(gcevec);
figure
handles(1) = loglog(lamvec,dirvec,'linewidth',3);
hold on
handles(2) = loglog(lamvec,gcdvec,'--','linewidth',3);
handles(3) = loglog(lamvec,eigvec,'r','linewidth',3);
handles(4) = loglog(lamvec,gcevec,'--r','linewidth',3);
loglog(lamvec(id),dirvec(id),'x','linewidth',3,'markersize',12)
loglog(lamvec(ig),gcdvec(ig),'x','linewidth',3,'markersize',12)
loglog(lamvec(ie),eigvec(ie),'r+','linewidth',3,'markersize',12)
loglog(lamvec(ic),gcevec(ic),'r+','linewidth',3,'markersize',12)
handles(5) = loglog(lamvec,gqrvec,'ok','linewidth',1);
title(sprintf('ep=%g,N=%d,M=%d',ep,N,M))
xlabel('\mu')
legend(handles,'Standard Basis','Standard GCV',...
       'Eigenfunction Basis','Eigenfunction GCV',...
       'Stable Basis','location','northwest')
hold off

% This can compute the lam=0 error, which I guess should be mu now that
% I think about it
GQR_reg.coef = Phi_reg\y;
yp = gqr_eval(GQR_reg,xx);
lam0_err = errcompute(yp,yy);

% Plot the results
figure
plot(x,y,'or')
hold on
plot(xx,yy,'linewidth',2)
plot(xx,eig_y_best,'--k','linewidth',2)
plot(xx,dir_y_best,'-.m','linewidth',2)
hold off
ylim([-1,2])
title(sprintf('ep=%g,lam=%g,N=%d,M=%d',ep,eig_lam_best,N,M))
legend('Data','True',...
    sprintf('Eig Opt err=%2.2g',eig_err_best),...
    sprintf('Stand Opt err=%2.2g',dir_err_best),...
    'location','south')
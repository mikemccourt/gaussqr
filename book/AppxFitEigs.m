% AppxFitEigs
% This script studies the quality of interpolation for increasingly many
% eigenfunctions used to approximate a low rank kernel matrix
% We study the Gaussian here, both because the rapid decay of its
% eigenvalus suggests a very low rank interpolation matrix and because we
% have a closed form for its eigenfunctions

% Define the Gaussian
rbf = @(e,r) exp(-(e*r).^2);
ep = .2;

% Define the data of interest
N = 100;
Neval = 500;
x = pickpoints(-1,1,N,'halton');
xeval = pickpoints(-1,1,Neval);
yf = @(x) cos(2*pi*x.^2);
y = yf(x);
yeval = yf(xeval);

% Form the full kernel matrix and determine its rank
K = rbf(ep,DistanceMatrix(x,x));
Keval = rbf(ep,DistanceMatrix(xeval,x));
rankK = rank(K);

% Precompute the Phi matrices
alpha = 1;
GQR = gqr_solveprep(1,x,ep,alpha,rankK+20);
Phi = gqr_phi(GQR,x);
Mvec = GQR.Marr;
lamvec = GQR.eig(Mvec);
rankvec = arrayfun(@(M)rank(bsxfun(@times,lamvec(1:M),Phi(:,1:M))*Phi(:,1:M)'),Mvec);
rankPhivec = arrayfun(@(M)rank(Phi(:,1:M)),Mvec);

plot(Mvec,rankvec)
hold on
plot(Mvec,rankK*ones(size(Mvec)),':')
plot(Mvec,rankPhivec,'--')
hold off
xlabel('M (number of eigenfunctions)')
ylabel('matrix rank')
legend('series kernel','closed form kernel','eigenfunctions','location','southeast')

% Loop through the eigenfunctions from 1 until slightly greater than the
% rank of the kernel matrix and follow the quality
errvec = arrayfun(@(M)errcompute(gqr_eval(gqr_rsolve(x,y,ep,alpha,M),xeval),yeval),Mvec);

% Compute the ranks at each step along the way
rankvec = arrayfun(@(M)rank(gqr_phi(1:M,x,ep,alpha)),Mvec);

% Plot the results
plotyy(Mvec,errvec,Mvec,rankvec,'loglog','semilogx')
% AppxFitEigs
% This script studies the quality of interpolation for increasingly many
% eigenfunctions used to approximate a low rank kernel matrix
% We study the Gaussian here, both because the rapid decay of its
% eigenvalus suggests a very low rank interpolation matrix and because we
% have a closed form for its eigenfunctions
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

% Define the Gaussian
rbf = @(e,r) exp(-(e*r).^2);
ep = .2;

% Define the data of interest
N = 100;
Neval = 500;
x = pickpoints(-1,1,N,'halton');
xeval = pickpoints(-1,1,Neval);
yf = @(x) 2*x.^3 + cos(2*pi*x.^2);
y = yf(x);
yeval = yf(xeval);

% Form the full kernel matrix and determine its rank
K = rbf(ep,DistanceMatrix(x,x));
Keval = rbf(ep,DistanceMatrix(xeval,x));
rankK = rank(K);

% Precompute the Phi matrices
alpha = 1;
GQR = gqr_solveprep(1,x,ep,alpha,rankK+20);
Mvec = GQR.Marr;
lamvec = GQR.eig(Mvec);
Phi = gqr_phi(GQR,x);
Phieval = gqr_phi(GQR,xeval);

% Create a function to evaluate the K series approximation with M terms
MtermK = @(M) bsxfun(@times,lamvec(1:M),Phi(:,1:M))*Phi(:,1:M)';
MtermKeval = @(M) bsxfun(@times,lamvec(1:M),Phieval(:,1:M))*Phi(:,1:M)';

% Study the rank of the Phi*Lambda*Phi' matrix and also the Phi matrix
rankvec = arrayfun(@(M)rank(MtermK(M)),Mvec);
rankPhivec = arrayfun(@(M)rank(Phi(:,1:M)),Mvec);

% Plot the results of the rank study
h_rank = figure;
plot(Mvec,rankvec)
hold on
plot(Mvec,rankK*ones(size(Mvec)),':')
plot(Mvec,rankPhivec,'--')
hold off
xlabel('M (number of eigenfunctions)')
ylabel('matrix rank')
legend('series kernel','closed form kernel','eigenfunctions','location','southeast')

% Compute the error associated with the K based interpolant
% Because K is a low rank matrix, we must use the psuedoinverse to produce
% a set of interpolation coefficients
errvec = arrayfun(@(M)errcompute(MtermKeval(M)*(pinv(MtermK(M))*y),yeval),Mvec);

% Loop through an increasing number of eigenfunctions and study the quality
% of the approximation
% This can be computed with the gaussqr toolbox functions:
%    errvec = arrayfun(@(M)errcompute(gqr_eval(gqr_rsolve(x,y,ep,alpha,M),xeval),yeval),Mvec);
% Below, we explicitly demonstrate the matrix computations
% Note that the least squares solution is obtained with K\y for an
% overdetermined system - this is the same notation for matrix inverse
errPhivec = arrayfun(@(M)errcompute(Phieval(:,1:M)*(Phi(:,1:M)\y),yeval),Mvec);

% Plot the results
h_err = figure;
semilogy(Mvec,errPhivec,'linewidth',2)
hold on
semilogy(Mvec,errvec,'--','linewidth',2)
hold off
xlabel('M - eigenfunctions used')
ylabel('2-norm of error')
legend('eigenfunction basis','low rank kernel','location','southwest')
% HSSVDPsivsK
% This example demonstrates the computation of the Psi matrix which appears
% in the HS-SVD and shows how the psi functions are related to the k
% functions through the eigenfunctions phi and a correction

% Define the points under consideration
N = 100;
x = pickpoints(-1,1,N,'cheb');

% Define the closed form of the Cinf kernel
K = @(x,z) .6 + .56*(.273 -.6*bsxfun(@plus,x.^2,z'.^2) + 1.27*x*z')./...
    (.8281 + .36*bsxfun(@plus,x.^2,z'.^2) - 1.308*x*z');

% Define the phi functions
% phi accepts row vector n and column vector x
%     it returns a length(x)-by-length(n) matrix
phi = @(n,x) bsxfun(@times,sqrt(2-(n==0)),cos(bsxfun(@times,n,acos(x))));
                
% Define the eigenvalues in two different ways
% We are fixing a = .4, b = .3
lam = @(n) (n==0)*.6 + (n>0).*(.3).^(n-1)*(.4*.7);

% Choose the number of eigenfunctions required to compute the HS-SVD
% Chosen arbitrarily here, but analysis could be used
M = N + 20;

% Evaluate the Phi matrix and diagonal of Lambda
% Again, the off-by-one issue pops up because the eigenvalues are indexed
% from 0 but Matlab starts indexing at 1
narr = 0:M-1;
Phi = phi(narr,x);
lamvec = lam(narr);

% Compute the closed form kernel matrix and its rank
Kcf = K(x,x);
rankK = rank(Kcf);

% Form the components of the HS-SVD
Phi1 = Phi(:,1:N);
Phi2 = Phi(:,N+1:end);
lamvec1 = lamvec(1:N);
lamvec2 = lamvec(N+1:end);

% Compute the Rbar matrix and the correction term
% Remember that the Lambda2 and inv(Lambda1) occur simultaneously
Rbar = bsxfun(@rdivide,lamvec2',lamvec1).*(Phi2'/Phi1');

% Form the full Psi matrix
Psi = Phi*[eye(N);Rbar];
rankPsi = rank(Psi);
Khssvd = bsxfun(@times,Psi,lamvec1)*Phi1';
err = norm(Kcf - Khssvd);

fprintf('rank(K) = %d, rank(Psi) = %d, error in HS-SVD = %g\n',rankK,rankPsi,err)

% Choose some points at which to plot the components of the HS-SVD
xeval = pickpoints(-1,1,1000);
Phieval = phi(narr,xeval);
Phieval1 = Phieval(:,1:N);
Phieval2 = Phieval(:,N+1:end);
Correction = Phieval2*Rbar;
Psieval = Phieval1 + Correction;

% Plot the eigenfunctions, correction and HS-SVD basis
h1 = figure;
plot(xeval,Phieval(:,97),':')
hold on
plot(xeval,Correction(:,97))
plot(xeval,Psieval(:,97),'--')
hold off
xlim([-.2 .2])
ylim([-2 2])
legend('eigenfunction','correction','hssvd basis')

h2 = figure;
plot(xeval,Phieval(:,98),':')
hold on
plot(xeval,Correction(:,98))
plot(xeval,Psieval(:,98),'--')
hold off
xlim([-.2 .2])
ylim([-2 2])
legend('eigenfunction','correction','hssvd basis')

h3 = figure;
plot(xeval,Phieval(:,99),':')
hold on
plot(xeval,Correction(:,99))
plot(xeval,Psieval(:,99),'--')
hold off
xlim([-.2 .2])
ylim([-2 2])
legend('eigenfunction','correction','hssvd basis')
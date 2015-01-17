% KernelsLowRank
% This example studies the rank of a kernel interpolation matrix as a
% function of the truncation length of the series used to generate it
% We consider only the Cinf Chebyshev kernel here

% Define the points under consideration
N = 100;
x = pickpoints(-1,1,N);

% Define the closed form of the Cinf kernel
K = @(x,z) .6 + .56*(.273 -.6*bsxfun(@plus,x.^2,z'.^2) + 1.27*x*z')./...
    (.8281 + .36*bsxfun(@plus,x.^2,z'.^2) - 1.308*x*z');

% Define the phi functions
% phi accepts row vector n and column vector x
%     it returns a length(x)-by-length(n) matrix
phi = @(n,x) bsxfun(@times,sqrt(2-(n==0)),cos(bsxfun(@times,n,acos(x))));
                
% Define the eigenvalues in two different ways
% We are fixing a = .4, b = .3, beta = 1
% By fixing beta = 1 we simplify the zeta function
%  computation to zeta(2) = pi^2/6
lam = @(n) (n==0)*.6 + (n>0).*(.3).^(n-1)*(.4*.7);

% Evaluate the Phix matrix and diagonal of Lambda
narr = 0:N-1;
Phi = phi(narr,x);
lamvec = lam(narr);

% Check the rank of the series interpolation matrix as a function of M
Mcheck = [1:49,50:50:N];
Philam = bsxfun(@times,lamvec,Phi);
rankvec = arrayfun(@(n) rank(Phi(:,1:n)*Philam(:,1:n)'),Mcheck);

% Compute the rank of the closed form kernel matrix
Kcf = K(x,x);

% Plot the results
semilogx(Mcheck,rankvec,'linewidth',3)
hold on
semilogx(Mcheck,ones(size(Mcheck))*rank(Kcf),'--','linewidth',2)
hold off
xlabel('series truncation length')
ylabel('rank')
legend('series kernel','closed form','location','southeast')
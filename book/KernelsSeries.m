% KernelsSeries
% This is an example of a kernel evaluated through series rather than
% through a closed form of the kernel
% The series is K(x,z) = sum lam_n phi_n(x) phi_n(z)
% The kernels of interest are the Chebyshev kernels of section 3.9.2
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.NORM_TYPE = 2;
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;

% Define the points under consideration
N = 10;
Neval = 500;
x = pickpoints(-1,1,Neval);
z = pickpoints(-1,1,N);

% Define the phi functions
% phi accepts row vector n and column vector x
%     it returns a length(x)-by-length(n) matrix
phi = @(n,x) bsxfun(@times,sqrt(2-(n==0)),...
                    cos(bsxfun(@times,n,acos(x))));
                
% Define the eigenvalues in two different ways
% We are fixing a = .4, b = .3, beta = 1
% By fixing beta = 1 we simplify the zeta function
%  computation to zeta(2) = pi^2/6
lamG = @(n) (n==0)*.6 + (n>0).*(.3).^(n-1)*(.4*.7);
lamA = @(n) (n==0)*.6 + (n>0).*.4./(n.^2*pi^2/6 + eps);

% Define the closed form of the kernels
% It will help us to define the degree 2 Bernoulli polynomial
% as well as an auxiliary term for the algebraic kernel
% Note these accept x and z as columns
KG = @(x,z) .6 + .56*(.273 -.6*bsxfun(@plus,x.^2,z'.^2) + 1.27*x*z')./...
    (.8281 + .36*bsxfun(@plus,x.^2,z'.^2) - 1.308*x*z');
B2 = @(x) x.^2 - x + 1/6;
aux = @(op,x,z) B2(abs(bsxfun(op,acos(x),acos(z)))/(2*pi));
KA = @(x,z) .6 + 6*(.4)*(aux(@plus,x,z') + aux(@minus,x,z'));

% Evaluate the close form kernels at the required data locations
KGcf = KG(x,z);
KAcf = KA(x,z);

% Evaluate the Phix and Phiz matrices and diagonal of Lambda
% Consider the first M=1000 terms, starting from zero
M = 1000;
narr = 0:M-1;
Phix = phi(narr,x);
Phiz = phi(narr,z);
lamvecG = lamG(narr);
lamvecA = lamA(narr);


% Compute the M-length series approximation to the kernels
% This is the Phix*Lambda*Phiz' computation
KGser = bsxfun(@times,Phix,lamvecG)*Phiz';
KAser = bsxfun(@times,Phix,lamvecA)*Phiz';

% Study the convergence of the series representation with increasing terms
% ncheck are the series lengths for which to check the quality
% Apply the lambda values to the kernel centers first to make the code
% simpler during the multiple evaluations for different M
% errcompute is a gaussqr function that returns how off the series is
ncheck = [1:49,50:50:M];
PhizlamG = bsxfun(@times,lamvecG,Phiz);
PhizlamA = bsxfun(@times,lamvecA,Phiz);
errvecG = arrayfun(@(n) errcompute(Phix(:,1:n)*PhizlamG(:,1:n)',KGcf),ncheck);
errvecA = arrayfun(@(n) errcompute(Phix(:,1:n)*PhizlamA(:,1:n)',KAcf),ncheck);

% Plot the error results
h_error = figure;
loglog(ncheck,errvecA,'linewidth',3)
hold on
loglog(ncheck,errvecG,'--','linewidth',3)
hold off
xlabel('summation length')
ylabel('2-norm series error')
legend(sprintf('\\beta=1 algebraic decay'),'geometric decay','location','northeast')

% Plot the kernels
h_kernels = figure;
plot(x,KAcf)
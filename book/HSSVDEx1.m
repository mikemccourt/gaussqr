% HSSVDEx1
% This example deals with Chebyshev kernels, and the use of the
% Hilbert-Schmidt SVD to stably interpolate with them
% These kernels were written up in a note by Fred Hickernell that I cannot
% find right now for the life of me
% As soon as I find it I will type up the explanation of it here

% Pick a function to sample from
% In this case we pick a mean function to test ourselves a bit
yf = @(x) 1./(1+25*x.^2);
% yf = @(x) 4*sin(3*x)./exp(3*x/2);

% Set up some data points at which to sample
N = 81;
x = pickpoints(-1,1,N,'cheb');
y = yf(x);

% Set up some evaluation points
NN = 100;
xx = pickpoints(-1,1,NN);
yy = yf(xx);

% Choose shape parameters for the kernel
% Not yet sure what the significance of these is
% We also may not have a closed form for the kernel so we will probably end
% up dealing with evaluation using the eigenfunction series
alpha = .5;
epvec = [logspace(-1,-.18,15),logspace(-.17,-.01,80)];

% Define the eigenfunctions and eigenvalues
% Note that lambdafunc takes in a vector and returns a vector
% phifunc takes in two vectors and returns a matrix
%    for size(a) = [1 p], size(b) = [q 1]
%    size(phifunc(a,b)) = [q p]
phifunc = @(n,x) sqrt(2)*cos(acos(x)*n);
lambdafunc = @(e,n) (n==0).*(1-alpha) + (n>0).*(alpha*(1-e)/e*e.^n);

% Determine how many eigenfunctions will be needed for accuracy
% The correct equation for this is
%                 M = N + ceil(log(eps)/log(ep));
% Unfortunately, this depends on ep, and I don't feel like reevaluating Phi
% everytime, so I'm going to figure out the max that will be needed and
% take them from the matrix Phi as needed within the loop below
Mmax = N + ceil(log(eps)/log(max(epvec)));

% Evaluate phi for all necessary eigenfunctions because they are unaffected
% by changes in epsilon
% Also note that the Phi1 matrix will never change, so we compute its
% inverse here, which I think is the same cost as factoring into L and U
% and reusing those factors down below a bunch of times
% We define nmax as varying from 0 to Mmax-1 so it still has M terms
nmax = 0:Mmax-1;
Phi = phifunc(nmax,x);
Phieval = phifunc(nmax,xx);
Phi1 = Phi(:,1:N);
Phieval1 = Phieval(:,1:N);

% Loop over the ep values under consideration
errvec = zeros(size(epvec));
errvecqr = zeros(size(epvec));
k = 1;
for ep=epvec
    % Evaluate the eigenvalues which are a function of epsilon
    M = N + ceil(log(eps)/log(max(ep)));
    Lambdavec = lambdafunc(ep,0:M-1);
    
    % First solve in the standard basis
    K = bsxfun(@times,Phi(:,1:M),Lambdavec)*Phi(:,1:M)';
    Keval = bsxfun(@times,Phieval(:,1:M),Lambdavec)*Phi(:,1:M)';
    warning('off','MATLAB:nearlySingularMatrix')
    yeval = Keval*(K\y);
    warning('on','MATLAB:nearlySingularMatrix')
    errvec(k) = errcompute(yeval,yy);

    % Now solve with the RBF-QR technique
    % Separate the Phi matrix into its blocks
    Phi2 = Phi(:,N+1:M);
    Phieval2 = Phieval(:,N+1:M);
    D = bsxfun(@rdivide,Lambdavec(N+1:M)',Lambdavec(1:N));
    Rbar = D.*(Phi2'/Phi1');
    Psi = Phi1 + Phi2*Rbar;
    Psieval = Phieval1 + Phieval2*Rbar;
    yevalqr = Psieval*(Psi\y);
    errvecqr(k) = errcompute(yevalqr,yy);
    
    % Store the interpolations for the first epsilon considered
    if k==1
        yplot = yeval;
        yplotqr = yevalqr;
    end
    
    k = k + 1;
end

% Plot the error as a function of epsilon
h_ep = figure;
loglog(epvec,errvecqr,'linewidth',3)
hold on
loglog(epvec,errvec,'r','linewidth',2)
hold off
legend('ChebQR','Standard basis')

% Plot the computed solutions
h_fit = figure;
plot(xx,yplot,'linewidth',2)
hold on
plot(xx,yplotqr,'g',xx,yy,':r','linewidth',3)
hold off
legend('Standard basis','ChebQR','True solution')
title(sprintf('N=%d, \\epsilon=%g',N,epvec(1)))

% Plot the pointwise error
h_err = figure;
err_ptwise = abs(yplot-yy);
errqr_ptwise = abs(yplotqr-yy)+eps;    % There seem to be a few places with 
                                        % zero error, so to make this look 
                                        % nicer a little shift by eps 
semilogy(xx,err_ptwise,'linewidth',2)
hold on
semilogy(xx,errqr_ptwise,'g','linewidth',3)
hold off
legend('Standard basis','ChebQR')
title(sprintf('Direct error = %g, QR error = %g',err,errqr))
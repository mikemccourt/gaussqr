% ParameterPower
% This example tests different methods for computing the power function of
% a Gaussian interpolation problem.  We use the Chebyshev kernel because we
% have access to its Mercer series

% Define the Chebyshev kernel with a = .5
Kf = @(b,x,z) .5 + (1-b)* ...
         (b*(1-b^2) - 2*b*bsxfun(@plus,x.^2,z.^2') + (1+3*b^2)*x*z')./ ...
         ((1-b^2)^2 + 4*b*(b*bsxfun(@plus,x.^2,z.^2')-(1+b^2)*x*z'));

% Define the eigenfunctions and eigenvalues
phifunc = @(n,x) sqrt(2)*cos(acos(x)*n);
lambdafunc = @(b,n) (n==0)*.5 + (n>0).*(.5*(1-b)/b*b.^n);

% Create some kernel centers
N = 11;
x = pickpoints(-1,1,N,'cheb');
Neval = 13;
xeval = pickpoints(-1,1,Neval,'halt');

% Define a region of shape parameters to test
bvec = logspace(-3,-.1,10);
nummat = zeros(length(Neval),length(bvec));
denvec = zeros(size(bvec));

warning('off','MATLAB:nearlySingularMatrix')

% Loop through and compute the power function
powvec = zeros(size(bvec));
powvecH = zeros(size(bvec));
powvecD = zeros(size(bvec));
k = 1;
for b=bvec
    % Try computing with the standard basis
    K = Kf(b,x,x);
    Keval = Kf(b,xeval,x);
    Kcenters = arrayfun(@(x)Kf(b,x,x),xeval);
    pvals = Kcenters - sum((Keval/K).*Keval,2);
    powvec(k) = sum(sqrt(abs(pvals)));
    
    % Try computing with the HS-SVD basis    
    M = N + ceil(log(eps)/log(b));
    n = 0:M-1;
    Phi = phifunc(n,x);
    Phieval = phifunc(n,xeval);
    Phi1 = Phi(:,1:N);
    Phi2 = Phi(:,N+1:end);
    lamvec = lambdafunc(b,n);
    lamvec1 = lamvec(1:N);
    lamvec2 = lamvec(N+1:end);
    CbarT = bsxfun(@rdivide,lamvec2',lamvec1).*(Phi2'/Phi1');
    Psi = Phi1 + Phi2*CbarT;
    Psieval = Phieval*[eye(N);CbarT];
    pvals = Kcenters - sum((Psieval/Psi).*Keval,2);
    powvecH(k) = sum(sqrt(abs(pvals)));
    
    % Try using the HS-SVD to compute the determinants
    % First the denominator
    logdetDenom = log(abs(det(Phi1))) + sum(log(lamvec1)) + log(abs(det(Psi)));
    denvec(k) = logdetDenom;
    % Then the numerator
    for m=1:Neval
        xp = [x;xeval(m)];
        Np = N + 1;
        M = Np + ceil(log(eps)/log(b));
        n = 0:M-1;
        Phi = phifunc(n,xp);
        Phi1 = Phi(:,1:Np);
        Phi2 = Phi(:,Np+1:end);
        lamvec = lambdafunc(b,n);
        lamvec1 = lamvec(1:Np);
        lamvec2 = lamvec(Np+1:end);
        CbarT = bsxfun(@rdivide,lamvec2',lamvec1).*(Phi2'/Phi1');
        Psi = Phi1 + Phi2*CbarT;
        logdetNum = log(abs(det(Phi1))) + sum(log(lamvec1)) + log(abs(det(Psi)));
        nummat(m,k) = logdetNum;
    end
    pvals = nummat(:,k) - denvec(k);
    powvecD(k) = sum(sqrt(exp(pvals)));
    
    k = k + 1;
end

warning('on','MATLAB:nearlySingularMatrix')

h_pow = figure;
h(3) = loglog(bvec,powvecD,'linewidth',2);
hold on
h(2) = loglog(bvec,powvecH,'--','linewidth',2);
h(1) = loglog(bvec,powvec,':','linewidth',2);
hold off
xlabel('b')
ylabel('power function 1-norm')
legend(h,'standard basis','HS-SVD basis','determinant','location','southeast')
ylim([1e-15 1])
% ex5c.m
% Comparing various methods of 2D BVP solvers
%  Trefethen's cheb matrix
%  Kansa's nonsymmetric
%  GaussQR regression

% The problem we are interested in solving is the Helmholtz equation:
%    Lap(u)+k^2 u = f
%    (x,y)=[-1,1]^2     homogeneous Dirichlet BC
%    f = exp(-10*((yy-1).^2+(xx-.5).^2));

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
GAUSSQR_PARAMETERS.ERROR_STYLE = 2; % Use absolute error

% f = @(x,y) exp(-10*((y-1).^2+(x-.5).^2));
% Solution below is u = (1-y^2)sin(pi*x)cosh(x+y)
fsol = @(x,y) (1-y.^2).*sin(pi*x).*cosh(x+y);
f = @(x,y) (1-y.^2).*(sin(pi*x).*cosh(x+y)*(1-pi^2) + 2*pi*cos(pi*x).*sinh(x+y)) - sin(pi*x).*((1+y.^2).*cosh(x+y)+4*y.*sinh(x+y)) + k^2*fsol(x,y);
% fsol = @(x,y) (1-y.^2).*sin(pi*x).*(1-exp(-10*((x-.5).^2+(y-.5).^2)));
% f = @(x,y) (y.^2-1).*sin(pi*x)*pi^2-2*sin(pi*x)+...
%     exp(-10*((x-.5).^2+(y-.5).^2)).*...
%     ((2*y.^2-2).*cos(pi*x)*pi.*(10-20*x)+...
%      40*(1-y.^2).*sin(pi*x)+...
%      (y.^2-1).*sin(pi*x).*(10-20*x).^2+...
%      (1-y.^2).*sin(pi*x)*pi^2+2*sin(pi*x)+...
%      4*y.*sin(pi*x).*(10-20*y)+...
%      (y.^2-1).*sin(pi*x).*(10-20*y).^2)
%     + k^2*fsol(x,y);

rbf = @(e,r) exp(-(e*r).^2);
drbf = @(e,r,dx) -2*e^2*dx.*exp(-(e*r).^2);
d2rbf = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);

% These are the functions needed for the Laplacian
Lrbf = @(e,r) 4*e^2*((e*r).^2-1).*exp(-(e*r).^2);

epvec = logspace(-1,1,20);

% Trefethen method first
N = 24;
[D,x] = cheb(N);
y = x;
[xx,yy] = meshgrid(x,y);
xx = xx(:); yy = yy(:);
usol = fsol(xx,yy);
b = find(abs(xx)==1 | abs(yy)==1); % Identify boundaries

D2 = D^2;
I = eye(N+1);
k = 9; % Helmholtz parameter
L = kron(I,D2) + kron(D2,I)+k^2*kron(I,I);
rhs = f(xx,yy);

L(b,:) = zeros(4*N,(N+1)^2); L(b,b) = eye(4*N);
rhs(b) = zeros(4*N,1);
u = L\rhs;
err_Trefethen = errcompute(u,usol);

% Kansa unsymmetric collocation with 1D kronecker products
r = DistanceMatrix(x,x);
errvec1D = [];
m = 1;
for ep=epvec
    A = rbf(ep,r);
    D2A = d2rbf(ep,r);
    D2 = D2A/A;
    L = kron(I,D2) + kron(D2,I) + k^2*kron(I,I);
    L(b,:) = zeros(4*N,(N+1)^2); L(b,b) = eye(4*N);
    u = L\rhs;
    errvec1D(m) = errcompute(u,usol);
    m = m+1;
end

pts = [xx,yy];
r = DistanceMatrix(pts,pts);
errvec2D = [];
m = 1;
for ep=epvec
    A = rbf(ep,r);
%     L = Lrbf(ep,r) + k^2*kron(I,I);
%     L(b,:) = A(b,:);
%     u = A*(L\rhs);
    LA = Lrbf(ep,r);
    L = LA/A + k^2*kron(I,I);
    L(b,:) = zeros(4*N,(N+1)^2); L(b,b) = eye(4*N);
    errvec2D(m) = errcompute(L\rhs,usol);
    m = m+1;
end

N = N+1;
m = 1;
errvecQ1D = [];
for ep=epvec
    [ep,alpha,Marr,lam] = rbfsolveprep(0,x,ep);
    phiMat = rbfphi(Marr,x,ep,alpha);
    phiMat2d = rbfphi(Marr,x,ep,alpha,2);
    [Q,R] = qr(phiMat);
    R1 = R(:,1:N);
    R2 = R(:,N+1:end);
    iRdiag = diag(1./diag(R1));
    R1s = iRdiag*R1;
    opts.UT = true;
    Rhat = linsolve(R1s,iRdiag*R2,opts);
    Ml = size(Marr,2);
    D = lam.^(repmat(sum(Marr(:,N+1:end),1)',1,N)-repmat(sum(Marr(:,1:N),1),Ml-N,1));
    Rbar = D.*Rhat';
    D2 = phiMat2d*[I;Rbar]/(phiMat*[I;Rbar]);
    L = kron(I,D2) + kron(D2,I) + k^2*kron(I,I);
    L(b,:) = zeros(4*(N-1),N^2); L(b,b) = eye(4*(N-1));
    errvecQ1D(m) = errcompute(L\rhs,usol);
    m = m+1;
end

GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .9;
m = 1;
errvecR1D = [];
for ep=epvec
    [ep,alpha,Marr] = rbfsolveprep(1,x,ep);
    phiMat = rbfphi(Marr,x,ep,alpha);
    phiMat2d = rbfphi(Marr,x,ep,alpha,2);
    D2 = phiMat2d/phiMat;
    L = kron(I,D2) + kron(D2,I) + k^2*kron(I,I);
    L(b,:) = zeros(4*(N-1),N^2); L(b,b) = eye(4*(N-1));
    errvecR1D(m) = errcompute(L\rhs,usol);
    m = m+1;
end

GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .5;
m = 1;
errvecR2D = [];
for ep=epvec
    [ep,alpha,Marr] = rbfsolveprep(1,pts,ep);
    phiMat = rbfphi(Marr,pts,ep,alpha);
    phiMat2d = rbfphi(Marr,pts,ep,alpha,[2,0])+rbfphi(Marr,pts,ep,alpha,[0,2]);
    L = phiMat2d/phiMat + k^2*kron(I,I);
    L(b,:) = zeros(4*(N-1),N^2); L(b,b) = eye(4*(N-1));
    errvecR2D(m) = errcompute(L\rhs,usol);
    m = m+1;
end

loglog(epvec,errvec1D,'--b','Linewidth',2),hold on
loglog(epvec,errvec2D,'--g','Linewidth',2)
loglog(epvec,errvecR1D,'b','Linewidth',3)
loglog(epvec,errvecR2D,'g','Linewidth',3)
loglog(epvec,errvecQ1D,'r','Linewidth',3)
loglog(epvec,err_Trefethen*ones(size(epvec)),'--k','Linewidth',2)
legend('kron collocation','2D collocation','kron regression','2D regression','kron QRsolve','Trefethen'),hold off




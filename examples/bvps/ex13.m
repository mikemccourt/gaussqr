% ex13.m
% Comparing various methods of 2D BVP solvers
% These methods all employ a tensor product structure
% differentiation matrix
%
% The methods we consider are
%  Trefethen's cheb matrix
%  Kansa's nonsymmetric
%  HS-SVD collocation
%
% The problem we are interested in solving is the Helmholtz equation:
%    Lap(u)+k^2 u = f
%    (x,y)=[-1,1]^2     Dirichlet BC
%
% There are a few different f choices we consider below

global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2; % Use absolute error
GAUSSQR_PARAMETERS.NORM_TYPE = inf; % Use absolute error

h_waitbar = waitbar(0,'Initializing');
warning('off','MATLAB:nearlySingularMatrix')

k = 7; % Helmholtz parameter

% f = @(x,y) exp(-10*((y-1).^2+(x-.5).^2));
% Solution below is u = (1-y^2)sin(pi*x)cosh(x+y)
% fsol = @(x,y) (1-y.^2).*sin(pi*x).*cosh(x+y);
% f = @(x,y) cosh(x+y).*((2+k^2-pi^2)*(1-y.^2).*sin(pi*x)-2*sin(pi*x)) + ...
%            sinh(x+y).*(2*pi*(1-y.^2).*cos(pi*x)-4*y.*sin(pi*x));
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
% fsol = @(x,y) 1./(1+x.^2+y.^2);
% f = @(x,y) 8*(x.^2+y.^2).*fsol(x,y).^3-4*fsol(x,y).^2+k^2*fsol(x,y);
fsol = @(x,y) besselj(0,6*sqrt(x.^2+y.^2));
f = @(x,y) (k^2-36)*besselj(0,6*sqrt(x.^2+y.^2));

rbf = @(e,r) exp(-(e*r).^2);
drbf = @(e,r,dx) -2*e^2*dx.*exp(-(e*r).^2);
d2rbf = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);

% These are the functions needed for the Laplacian
Lrbf = @(e,r) 4*e^2*((e*r).^2-1).*exp(-(e*r).^2);

N = 20;
epvec = logspace(-1,1,35);

waitbar(0,h_waitbar,'Solving Trefethen System');

% Trefethen method first
[D,x] = cheb(N);
y = x;
[xx,yy] = meshgrid(x,y);
xx = xx(:); yy = yy(:);
usol = fsol(xx,yy);
b = (abs(xx)==1 | abs(yy)==1); % Identify boundaries

D2 = D^2;
I = eye(N);
L = kron(I,D2) + kron(D2,I) + k^2*kron(I,I);
rhs = f(xx,yy);

L(b,:) = zeros(4*(N-1),N^2); L(b,b) = eye(4*(N-1));
rhs(b) = zeros(4*(N-1),1);rhs(b) = usol(b);
u = L\rhs;
err_Trefethen = errcompute(u,usol);

% Kansa unsymmetric collocation with 1D kronecker products
r = DistanceMatrix(x,x);
errvec1D = zeros(size(epvec));
m = 1;
for ep=epvec
    waitbar(m/(2*length(epvec)),h_waitbar,sprintf('Solving Standard Basis, ep=%3.2f',ep));
    
    A = rbf(ep,r);
    D2A = d2rbf(ep,r);
    D2 = D2A/A;
    L = kron(I,D2) + kron(D2,I) + k^2*kron(I,I);
    L(b,:) = zeros(4*(N-1),N^2); L(b,b) = eye(4*(N-1));
    u = L\rhs;
    errvec1D(m) = errcompute(u,usol);
    
    m = m+1;
end

m = 1;
errvecQ1D = zeros(size(epvec));
alpha = 1;
for ep=epvec
    waitbar((m + length(epvec))/(2*length(epvec)),h_waitbar,sprintf('Solving HS-SVD Basis, ep=%3.2f',ep));
    
    GQR = gqr_solveprep(0,x,ep,alpha);
    phiMat_1 = GQR.stored_phi1;
    phiMat_2 = GQR.stored_phi2;
    phiMat2d = gqr_phi(GQR,x,2);
    Rbar = GQR.Rbar;
    D2 = phiMat2d*[I;Rbar]/(phiMat_1+phiMat_2*Rbar);
    L = kron(I,D2) + kron(D2,I) + k^2*kron(I,I);
    L(b,:) = zeros(4*(N-1),N^2); L(b,b) = eye(4*(N-1));
    errvecQ1D(m) = errcompute(L\rhs,usol);
    
    m = m+1;
end

waitbar(1,h_waitbar,'Plotting');

h = figure;
loglog(epvec,errvecQ1D,'r','Linewidth',3),hold on
loglog(epvec,errvec1D,'b','Linewidth',2)
loglog(epvec,err_Trefethen*ones(size(epvec)),'--k','Linewidth',2)
ylim([1e-13 1e1])
legend('HS-SVD','Fasshauer','Trefethen','location','west'),hold off
xlabel('\epsilon')
ylabel('Absolute sup-norm error')

close(h_waitbar)
warning('on','MATLAB:nearlySingularMatrix')
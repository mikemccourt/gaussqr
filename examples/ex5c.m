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
Mextramax = GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC;
Mfactor = GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC;

f = @(x,y) exp(-10*((y-1).^2+(x-.5).^2));
% Solution below is u = (1-y^2)sin(pi*x)cosh(x+y)
% f = @(x,y) (1-y.^2).*(sin(pi*x).*cosh(x+y)*(1-pi^2) + 2*pi*cos(pi*x).*sinh(x+y)) - sin(pi*x)*((1+y.^2).*cosh(x+y)+4*y.*sinh(x+y)) + k^2*(1-y.^2)*sin(pi*x).*cosh(x+y);
rbf = @(e,r) exp(-(e*r).^2);
d2rbf = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);
Lrbf = @(e,r) 4*e^2*((e*r).^2-1).*exp(-(e*r).^2);

epvec = logspace(-1,1,20);

% Trefethen method first
N = 24;
[D,x] = cheb(N);
y = x;
[xx,yy] = meshgrid(x,y);
xx = xx(:); yy = yy(:);
rhs = f(xx,yy);
D2 = D^2;
I = eye(N+1);
k = 9;
L = kron(I,D2) + kron(D2,I)+k^2*kron(I,I);
b = find(abs(xx)==1 | abs(yy)==1);
L(b,:) = zeros(4*N,(N+1)^2); L(b,b) = eye(4*N);
rhs(b) = zeros(4*N,1);
u = L\rhs;
uu = reshape(u,N+1,N+1);
val_Trefethen = uu(N/2+1,N/2+1);

% Kansa unsymmetric collocation
r = DistanceMatrix(x,x);
% for ep=epvec
%     A = rbf(ep,r);
%     D2A = d2rbf(ep,r);
% end
pts = [xx,yy];
r = DistanceMatrix(pts,pts);
Amat = rbf(1,r);
Lmat = Lrbf(1,r);
Dmat = Lmat/Amat;
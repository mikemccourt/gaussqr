% ex12.m
% This example shows the effect of ill-conditioning on derivatives
% This considers a few different methods of approximation
%
% Specifically, the function we are considering here is
%   u(x,y) = sinh(x)cosh(y)
% And we consider different methods of approximating the Laplacian
% on the domain [-1,1]^2
% The Laplacian of that function is Lap(u) = 2u
%
% The methods we consider are:
%   Trefethen: Cheb
%   Collocation 2D
%   Collocation 1D
%   GaussQR 2D
%   GaussQRr 2D

rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2; % Use absolute error

ufunc = @(x,y) sinh(x).*cosh(y);
Lf = @(x,y) 2*sinh(x).*cosh(y);
N = 24;

rbf = @(e,r) exp(-(e*r).^2);
drbf = @(e,r,dx) -2*e^2*dx.*exp(-(e*r).^2);
d2rbf = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);
Lrbf = @(e,r) 4*e^2*((e*r).^2-1).*exp(-(e*r).^2);

[D,x] = cheb(N);
y = x;
[xx,yy] = meshgrid(x,y);
xx = xx(:); yy = yy(:);
u = ufunc(xx,yy);
Lfu = ufunc(xx,yy);
D2 = D^2;
I = eye(N+1);
L = kron(I,D2) + kron(D2,I);
err_Trefethen = errcompute(L*u,Lfu);

epvec = logspace(-1,1,40);
pts = [xx,yy];
rp = DistanceMatrix(pts,pts);
r = DistanceMatrix(x,x);
errvec1D = [];
errvec2D = [];
k = 1;
for ep=epvec
  Amat = rbf(ep,rp);
  b = Amat\u;
  Lmat = Lrbf(ep,rp);
  errvec2D(k) = errcompute(Lmat*b,Lfu);
  A = rbf(ep,r);
  d2A = d2rbf(ep,r);
  I = eye(size(r));
  D = d2A/A;
  L = kron(I,D) + kron(D,I);
  errvec1D(k) = errcompute(L*u,Lfu);
  k = k + 1;
end

N = size(x,1);
epvec = logspace(-1,1,40);
errvecR2d = [];
% errvecQ2d = [];
errvecR1d = [];
errvecQ1d = [];
k = 1;
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
  L = kron(I,D2) + kron(D2,I);
  errvecQ1d(k) = errcompute(L*u,Lfu);
  
  [ep,alpha,Marr] = rbfsolveprep(1,x,ep);
  phiMat = rbfphi(Marr,x,ep,alpha);
  phiMat2d = rbfphi(Marr,x,ep,alpha,2);
  D2 = phiMat2d/phiMat;
  L = kron(I,D2) + kron(D2,I);
  errvecR1d(k) = errcompute(L*u,Lfu);
  
  % Note that alpha was defined earlier in rbfsolveprep
  GQR = rbfqrr_solve(pts,u,ep,alpha);
  Lu = rbfqr_eval(GQR,pts,[2,0]) + rbfqr_eval(GQR,pts,[0,2]);
  errvecR2d(k) = errcompute(Lu,Lfu);
  
%   GQR = rbfqr_solve(pts,u,ep,alpha);
%   Lu = rbfqr_eval(GQR,pts,[2,0]) + rbfqr_eval(GQR,pts,[0,2]);
%   errvecQ2d(k) = errcompute(Lu,Lfu);
  
  k = k + 1;
end

loglog(epvec,errvec1D,':b','LineWidth',2),hold on
loglog(epvec,errvec2D,':g','LineWidth',2)
loglog(epvec,errvecR1D,'--b','LineWidth',2)
loglog(epvec,errvecR2D,'--g','LineWidth',2)
loglog(epvec,errvecQ1D,'--r','LineWidth',2)
loglog(epvec,err_Trefethen*ones(size(epvec)),'--k','LineWidth',2),hold off
legend('kron Collocation','2D Collocation','kron Regression','kron QRsolve','2D Regression','Trefethen')

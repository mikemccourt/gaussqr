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
%   GaussQRr 2D

rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2; % Use absolute error

N = 24;
ufunc = @(x,y) sinh(x).*cosh(y);
Lf = @(x,y) 2*sinh(x).*cosh(y);
ufunc = @(x,y) sin(x+y)+exp(-((x-.5).^2+(y+.25).^2));
Lf = @(x,y) -2*sin(x+y)+(-4+(1-2*x).^2+(.5+2*y).^2).*exp(-((x-.5).^2+(y+.25).^2));
ufunc = @(x,y) 1./(1+x.^2+y.^2);
Lf = @(x,y) 4*(-1+x.^2+y.^2)./(1+x.^2+y.^2).^3;

rbf = @(e,r) exp(-(e*r).^2);
drbf = @(e,r,dx) -2*e^2*dx.*exp(-(e*r).^2);
d2rbf = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);
Lrbf = @(e,r) 4*e^2*((e*r).^2-1).*exp(-(e*r).^2);

[D,x] = cheb(N);
y = x;
[xx,yy] = meshgrid(x,y);
xx = xx(:); yy = yy(:);
u = ufunc(xx,yy);
Lfu = Lf(xx,yy);
D2 = D^2;
I = eye(N+1);
L = kron(I,D2) + kron(D2,I);
err_Trefethen = errcompute(L*u,Lfu);

epvec = logspace(-2,1,40);
pts = [xx,yy];
rp = DistanceMatrix(pts,pts);
r = DistanceMatrix(x,x);
errvec1d = [];
errvec2d = [];
k = 1;
for ep=epvec
  Amat = rbf(ep,rp);
  b = Amat\u;
  Lmat = Lrbf(ep,rp);
  errvec2d(k) = errcompute(Lmat*b,Lfu);
  A = rbf(ep,r);
  d2A = d2rbf(ep,r);
  I = eye(size(r));
  D = d2A/A;
  L = kron(I,D) + kron(D,I);
  errvec1d(k) = errcompute(L*u,Lfu);
  k = k + 1;
end

N = size(x,1);
errvecR2d = [];
errvecR1d = [];
errvecQ1d = [];
k = 1;
alpha = 3;
for ep=epvec
  [ep,alpha,Marr,lam] = gqr_solveprep(0,x,ep,alpha);
  phiMat = gqr_phi(Marr,x,ep,alpha);
  phiMat2d = gqr_phi(Marr,x,ep,alpha,2);
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
  
  [ep,alpha,Marr] = gqr_solveprep(1,x,ep,alpha);
  phiMat = gqr_phi(Marr,x,ep,alpha);
  phiMat2d = gqr_phi(Marr,x,ep,alpha,2);
  D2 = phiMat2d/phiMat;
  L = kron(I,D2) + kron(D2,I);
  errvecR1d(k) = errcompute(L*u,Lfu);
  
  % Note that alpha was defined earlier in gqr_solveprep
  GQR = gqr_rsolve(pts,u,ep,alpha);
  Lu = gqr_eval(GQR,pts,[2,0]) + gqr_eval(GQR,pts,[0,2]);
  errvecR2d(k) = errcompute(Lu,Lfu);
  
  k = k + 1;
end

loglog(epvec,errvecR1d,'b','LineWidth',3),hold on
loglog(epvec,errvecR2d,'g','LineWidth',3)
loglog(epvec,errvecQ1d,'r','LineWidth',3)
loglog(epvec,errvec1d,'-.b','LineWidth',2)
loglog(epvec,errvec2d,'-.g','LineWidth',2)
loglog(epvec,err_Trefethen*ones(size(epvec)),'--k','LineWidth',2),hold off
legend('kron Regression','2D Regression','kron QRsolve','kron Collocation',...
       '2D Collocation','Trefethen','Location','Northwest')

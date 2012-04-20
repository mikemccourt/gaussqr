% ex12b.m
% This example shows the order of convergence as a function of the
% derivative which is being approximated
%
% Specifically, the function we will consider here is
%   u(x) = x*sin(x)
% And we consider different orders of approximation on the domain [-3,3]^2
%
% We will be using primarily GaussQR, until the size of the problem demands
% the regression

rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2; % Use absolute error
GAUSSQR_PARAMETERS.NORM_TYPE = inf; % Use inf norm

Nvec = [6,8,10,12,14,16,18,20];
NvecR = 3*Nvec;
% Nvec = [10,15,20,25,30,35,40,45,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190];
NN = 200;
fopt = 4;
alpha = 2;

uf = @(x) x.*sin(x);
ufd = @(x) sin(x) + x.*cos(x);
ufd2 = @(x) 2*cos(x)-x.*sin(x);
ufd3 = @(x) -3*sin(x)-x.*cos(x);
ufd4 = @(x) x.*sin(x)-4*cos(x);

% uf = @(x) exp(sin(x));
% ufd = @(x) exp(sin(x)).*cos(x);
% ufd2 = @(x) exp(sin(x)).*(cos(x).^2-sin(x));
% ufd3 = @(x) exp(sin(x)).*cos(x).*(-3*sin(x)+cos(x).^2-1);
% ufd4 = @(x) exp(sin(x)).*(sin(x).*(3*sin(x)+1)+cos(x).^4-2*cos(x).^2.*(3*sin(x)+2));

% uf = @(x) 1./(1+x.^2);
% ufd = @(x) -2*x./(1+x.^2).^2;
% ufd2 = @(x) (6*x.^2-2)./(1+x.^2).^3;
% ufd3 = @(x) -(24*x.*(x.^2-1))./(1+x.^2).^4;
% ufd4 = @(x) 24*(5*x.^4-10*x.^2+1)./(1+x.^2).^5;

% Consider accuracy for a fixed kernel
ep = .1;

errvecQ = [];
errvecQ1d = [];
errvecQ2d = [];
errvecQ3d = [];
errvecQ4d = [];
k = 1;
for N=Nvec
    x = pickpoints(-3,3,N,'cheb');
    xx = pickpoints(-3,3,NN);
    rhs = uf(x);
    u = uf(xx);u1 = ufd(xx);u2 = ufd2(xx);u3 = ufd3(xx);u4 = ufd4(xx);
    I = eye(N);
    
    [ep,alpha,Marr,lam] = rbfsolveprep(0,x,ep,alpha);
    phiMat = rbfphi(Marr,x,ep,alpha);
    phiMat0d = rbfphi(Marr,xx,ep,alpha,0);
    phiMat1d = rbfphi(Marr,xx,ep,alpha,1);
    phiMat2d = rbfphi(Marr,xx,ep,alpha,2);
    phiMat3d = rbfphi(Marr,xx,ep,alpha,3);
    phiMat4d = rbfphi(Marr,xx,ep,alpha,4);
    
    [Q,R] = qr(phiMat);
    R1 = R(:,1:N);
    R2 = R(:,N+1:end);
    opts.UT = true;
    Rhat = linsolve(R1,R2,opts);
    Ml = size(Marr,2);
    L = lam.^(repmat(sum(Marr(:,N+1:end),1)',1,N)-repmat(sum(Marr(:,1:N),1),Ml-N,1));
    Rbar = L.*Rhat';
    b = (phiMat*[I;Rbar])\rhs;
    errvecQ(k) = errcompute(phiMat0d*([I;Rbar]*b),u);
    errvecQ1d(k) = errcompute(phiMat1d*([I;Rbar]*b),u1);
    errvecQ2d(k) = errcompute(phiMat2d*([I;Rbar]*b),u2);
    errvecQ3d(k) = errcompute(phiMat3d*([I;Rbar]*b),u3);
    errvecQ4d(k) = errcompute(phiMat4d*([I;Rbar]*b),u4);
  
    k = k + 1;
end

GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .3;

errvecR = [];
errvecR1d = [];
errvecR2d = [];
errvecR3d = [];
errvecR4d = [];
k = 1;
for N=NvecR
    x = pickpoints(-3,3,N,'cheb');
    xx = pickpoints(-3,3,NN);
    rhs = uf(x);
    u = uf(xx);u1 = ufd(xx);u2 = ufd2(xx);u3 = ufd3(xx);u4 = ufd4(xx);
    
    [ep,alpha,Marr] = rbfsolveprep(1,x,ep,alpha);
    phiMat = rbfphi(Marr,x,ep,alpha);
    phiMat0d = rbfphi(Marr,xx,ep,alpha,0);
    phiMat1d = rbfphi(Marr,xx,ep,alpha,1);
    phiMat2d = rbfphi(Marr,xx,ep,alpha,2);
    phiMat3d = rbfphi(Marr,xx,ep,alpha,3);
    phiMat4d = rbfphi(Marr,xx,ep,alpha,4);
    
    b = phiMat\rhs;
    errvecR(k) = errcompute(phiMat0d*b,u);
    errvecR1d(k) = errcompute(phiMat1d*b,u1);
    errvecR2d(k) = errcompute(phiMat2d*b,u2);
    errvecR3d(k) = errcompute(phiMat3d*b,u3);
    errvecR4d(k) = errcompute(phiMat4d*b,u4);
  
    k = k + 1;
end

semilogy(Nvec,[errvecQ;errvecQ1d;errvecQ2d;errvecQ3d;errvecQ4d]','linewidth',2)
ylim([1e-16 1e1])
xlim([min(Nvec),max(Nvec)])
xlabel('input points N')
ylabel('RMS relative error')
legend('interp','1 deriv','2 deriv','3 deriv','4 deriv')
title('Full series approximation')
figure
semilogy(NvecR,[errvecR;errvecR1d;errvecR2d;errvecR3d;errvecR4d]','linewidth',2)
ylim([1e-16 1e1])
xlim([min(NvecR),max(NvecR)])
xlabel('input points N')
ylabel('RMS relative error')
legend('interp','1 deriv','2 deriv','3 deriv','4 deriv')
title('Low-rank series approximation')

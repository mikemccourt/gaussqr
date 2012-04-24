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
GAUSSQR_PARAMETERS.ERROR_STYLE = 4; % Use absolute error

Nvec = [6,8,10,12,14,16,18,20,22];
NvecR = 3*Nvec;
% Nvec = [10,15,20,25,30,35,40,45,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190];
NN = 200;
xx = pickpoints(-3,3,NN);
fopt = 3;

switch fopt
    case 1
        uf = @(x,n) (n==0)*(x.*sin(x))+...
            (n==1)*(sin(x)+x.*cos(x))+...
            (n==2)*(2*cos(x)-x.*sin(x))+...
            (n==3)*(-3*sin(x)-x.*cos(x))+...
            (n==4)*(x.*sin(x)-4*cos(x));
    case 2
        uf = @(x,n) (n==0)*(exp(sin(x)))+...
            (n==1)*(exp(sin(x)).*cos(x))+...
            (n==2)*(exp(sin(x)).*(cos(x).^2-sin(x)))+...
            (n==3)*(exp(sin(x)).*cos(x).*(-3*sin(x)+cos(x).^2-1))+...
            (n==4)*(exp(sin(x)).*(sin(x).*(3*sin(x)+1)+cos(x).^4-2*cos(x).^2.*(3*sin(x)+2)));
    case 3
        uf = @(x,n) (n==0)*(1./(1+x.^2))+...
            (n==1)*(-2*x./(1+x.^2).^2)+...
            (n==2)*((6*x.^2-2)./(1+x.^2).^3)+...
            (n==3)*(-(24*x.*(x.^2-1))./(1+x.^2).^4)+...
            (n==4)*(24*(5*x.^4-10*x.^2+1)./(1+x.^2).^5);
end

% Consider accuracy for a fixed kernel
ep = .1;
alpha = 2;
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .4;

errmatQ = zeros(5,length(Nvec));
errmatR = zeros(5,length(Nvec));
k = 1;
for N=Nvec
    x = pickpoints(-3,3,N,'cheb');
    xR = pickpoints(-3,3,3*N,'cheb');
    
    GQR = rbfqr_solve(x,uf(x,0),ep,alpha);
    GQRr = rbfqrr_solve(xR,uf(xR,0),ep,alpha);
    for d=0:4
        errmatQ(d+1,k) = errcompute(rbfqr_eval(GQR,xx,d),uf(xx,d));
        errmatR(d+1,k) = errcompute(rbfqr_eval(GQRr,xx,d),uf(xx,d));
    end
  
    k = k + 1;
end

semilogy(Nvec,errmatQ,'linewidth',2)
ylim([1e-16 1e1])
xlim([min(Nvec),max(Nvec)])
xlabel('input points N')
ylabel('RMS relative error')
legend('interp','1 deriv','2 deriv','3 deriv','4 deriv')
title('Full series approximation')

figure
semilogy(3*Nvec,errmatR,'linewidth',2)
ylim([1e-16 1e1])
xlim([min(NvecR),max(NvecR)])
xlabel('input points N')
ylabel('RMS relative error')
legend('interp','1 deriv','2 deriv','3 deriv','4 deriv')
title('Low-rank series approximation')

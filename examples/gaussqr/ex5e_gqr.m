% ex5e
% This problem makes a comparison between the finite difference method and
% the spectral collocation method
% This is just to produce some pictures for a talk

rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = 1;

rhs_func = @(x) exp(x).*tan(x);

Nvec = 2.^(3:10)+1;

errvals = zeros(5,9);

x5 = pickpoints(0,1,5);
N2 = 16;
A5 = [1 0 0 0 0;N2 -2*N2+1 N2 0 0;0 N2 -2*N2+1 N2 0;0 0 N2 -2*N2+1 N2;0 0 0 0 1];
b5 = [1;exp(x5(2:end-1)).*tan(x5(2:end-1));0];
u5 = A5\b5;

errvals(:,1) = u5;

k = 2;
m = 1;
for N=Nvec
    x = pickpoints(0,1,N);
    A = (N-1)^2*gallery('tridiag',N)+eye(N);
    A([1,end],:) = zeros(size(A([1,end],:)));
    A([1,end],[1,end]) = eye(size(A([1,end],[1,end])));
    rhs = rhs_func(x);
    rhs([1,end]) = [1;0];
    
    u = A\rhs;
    i = 1;
    for j=1:5
        errvals(j,m) = u(i);
        i = i + k;
    end
    k = k*2;
    m = m + 1;
end

de = diff(errvals,1,2);
orderFD = log2(abs(de(2:end-1,1:end-1)./de(2:end-1,2:end)));
plot(x,u,'--b','linewidth',2),hold on,plot(x5,u5,'or','linewidth',3),hold off
ylim([0 1.01])
legend('True solution','N=5 FD solution','location','northeast')
xlabel('x')
ylabel('u')
title('u``+u=exp(x)tan(x), u(0)=1, u(1)=0')


% N=1000;
% x = 3*(2*rand(N,1)-1);
% y = 3*(2*rand(N,1)-1);
% s = find(abs(y)<=abs(x).^(-3) & abs(y)<=abs(x).^(-1) & abs(y)<=abs(x).^(-.33));
% xt = pickpoints(-3,3,floor(.1*N),'inner')+eps;
% yt = min([abs(xt).^(-3),abs(xt).^(-1),abs(xt).^(-.33),3*ones(size(xt))],[],2);
% x = [x(s);xt;xt];
% y = [y(s);yt;-yt];
% plot(x,y+3,'b.'),hold on
% plot(x+3,y,'b.')
% plot(x,y-3,'b.'),hold off
% 
% TRI = delaunay(x,y);
% k = 1;
% dump = [];
% while k<size(TRI,1)
%     triplot(TRI([1:k-1,k+1:end],:),x,y),hold on,triplot(TRI(k,:),x,y,'r'),hold off
%     keepit = input('good? ')
%     if keepit==2
%     elseif keepit==3
% %         TRI(k,:) = [];
%         dump = [dump,k];
%     else
%         break
%     end
%     k = k + 1;
% end

f = @(x) abs(x).*cos(x).^3;
Nvec = 10*2.^(0:5);
NN = 1000;
xx = pickpoints(-1,1,NN);
yy = f(xx);

errmat = zeros(6,5);
k = 1;
for N=Nvec
    x = pickpoints(-1,1,N);
    y = f(x);
    for m=1:5
        spl = spapi(m+1,x,y);
        yp = fnval(spl,xx);
        errmat(k,m) = errcompute(yp,yy);
    end
    k = k + 1;
end

f = @(x) exp(-x.^2);
Nvec = 4:2:30;
NN = 200;
xx = pickpoints(-1,1,NN);
yy = f(xx);
ep = .0001;
alpha = 1;

errvec = zeros(size(Nvec));
k = 1;
for N=Nvec
    x = pickpoints(-1,1,N,'cheb');
    y = f(x);
%     GQR = gqr_solve(x,y,ep,alpha);
%     yp = gqr_eval(GQR,xx);
    po = polyfit(x,y,N-1);
    yp = polyval(po,xx);
%     A = zeros(N);
%     AA = zeros(NN,N);
%     for m=1:N
%         A(:,m) = x.^(m-1);
%         AA(:,m) = xx.^(m-1);
%     end
%     c = A\y;
%     yp = AA*c;
    errvec(k) = errcompute(yp,yy);
    k = k + 1;
end

semilogy(Nvec,errvec,'linewidth',3)
ylim([1e-17 1e-2])
xlim([min(Nvec),max(Nvec)])
xlabel('N - data points')
ylabel('Absolute error')
title('Exponential convergence (until machine precision)')



fsol = @(x) cosh(x)+x;
frhs = @(x) cosh(x);
Nvec = 4:2:30;
NN = 200;
xx = pickpoints(0,1,NN);
yy = fsol(xx);
ep = .0001;
alpha = 1;

errvec = zeros(size(Nvec));
k = 1;
for N=Nvec
    x = pickpoints(0,1,N,'cheb');
    y = frhs(x);
    GQR = gqr_rsolve(x,y,ep,alpha);
    phi = gqr_phi(GQR.Marr,x,ep,alpha);
    phi2 = gqr_phi(GQR.Marr,x,ep,alpha,2);
    A = [phi(1,:);phi2(2:end-1,:);phi(end,:)];
    y = [fsol(x(1));frhs(x(2:end-1));fsol(x(end))];
    GQR.coef = A\y;
    yp = gqr_eval(GQR,xx);
    errvec(k) = errcompute(yp,yy);
    k = k + 1;
end

semilogy(Nvec,errvec,'linewidth',3)
ylim([1e-17 1e-2])
xlim([min(Nvec),max(Nvec)])
xlabel('N - data points')
ylabel('Absolute error')
title('u``=cosh(x), u(x)=cosh(x)+x')




fsol = @(x) 1./(1+x.^2);
frhs = @(x) -2./(1+x.^2).^2+8*x.^2./(1+x.^2).^3;
N = 35;
NN = 200;
x = pickpoints(-2,2,N,'cheb');
y = frhs(x);
xx = pickpoints(-2,2,NN);
yy = fsol(xx);
epvec = logspace(-1,.6,30);
alpha = 2;

DM_DATA = DistanceMatrix(x,x);
DM_EV = DistanceMatrix(xx,x);
rbf = @(ep,r) exp(-(ep*r).^2);
rbf2 = @(ep,r) (4*ep^4*r.^2-2*ep^2).*exp(-(ep*r).^2);

errvec = zeros(size(epvec));
errdir = zeros(size(epvec));
k = 1;
for ep=epvec
%     GQR = gqr_rsolve(x,y,ep,alpha);
%     phi = gqr_phi(GQR.Marr,x,ep,alpha);
%     phi2 = gqr_phi(GQR.Marr,x,ep,alpha,2);
%     A = [phi(1,:);phi2(2:end-1,:);phi(end,:)];
    GQR = gqr_solve(x,y,ep,alpha);
    phi = gqr_phi(GQR.Marr,x,ep,alpha);
    phi2 = gqr_phi(GQR.Marr,x,ep,alpha,2);
    A = [phi(1,:);phi2(2:end-1,:);phi(end,:)]*[eye(N);GQR.Rbar];
    rhs = [fsol(x(1));frhs(x(2:end-1));fsol(x(end))];
    GQR.coef = A\rhs;
    yp = gqr_eval(GQR,xx);
    errvec(k) = errcompute(yp,yy);
    
    IM = rbf(ep,DM_DATA);
    IM2 = rbf2(ep,DM_DATA);
    EM = rbf(ep,DM_EV);
    A = [IM(1,:);IM2(2:end-1,:);IM(end,:)];
    yp = EM *(A\rhs);
    errdir(k) = errcompute(yp,yy);
    
    k = k + 1;
end
    
D = cheb(N,-2,2);
D2 = D*D;
A = [[1,zeros(1,N-1)];D2(2:end-1,:);[zeros(1,N-1),1]];
c = A\rhs;
po = polyfit(x,c,N-1);
yp = polyval(po,xx);
errcheb = errcompute(yp,yy);

loglog(epvec,errvec,'linewidth',3),hold on
loglog(epvec,errdir,'-xg','linewidth',2)
loglog(epvec,errcheb*ones(size(epvec)),'--r','linewidth',2),hold off
%ylim([1e-17 1e-2])
xlim([min(epvec),max(epvec)])
xlabel('\epsilon')
ylabel('Absolute error')
title('u``=f, u(x)=1/(1+x^2)')
legend('True RBF','RBF-Direct','Polynomial','Location','Northeast')







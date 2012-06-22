% ex5b2.m
% Comparing various methods of 2-pt BVP solvers
%  Trefethen's cheb matrix
%  Kansa's nonsymmetric
%  GaussQR regression (M<N)

% The problem we are interested in solving is
%    u_xx = -9*pi^2*sin(3*pi*x)-pi^2*cos(pi*x)  u(-1)=u(1)=0

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
GAUSSQR_PARAMETERS.ERROR_STYLE = 2; % Use absolute error
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .5;

f = @(x) -9*pi^2*sin(3*pi*x)-pi^2*cos(pi*x);
exact = @(x) sin(3*pi*x)+cos(pi*x)+1;

rbf = @(e,r) exp(-(e*r).^2);
d2rbf = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);

N = 80;
NN = 200;
xx = pickpoints(-1,1,NN);
uu = exact(xx);

warning off MATLAB:nearlySingularMatrix

% Trefethen's method
[D,t] = cheb(N);
D2 = D^2;
D2([1,end],:) = [1,zeros(1,N-1);zeros(1,N-1),1];
rhs = [exact(t(1));f(t(2:end-1));exact(t(end))];
u = D2\rhs;
[pcoef,S,mu] = polyfit(t,u,N-1);
err_Trefethen = errcompute(polyval(pcoef,xx,S,mu),uu)

% Kansa's method
epvec = logspace(-2,2,31);
k = 1;
err_Kansa = [];
r = DistanceMatrix(t,t);
reval = DistanceMatrix(xx,t);
rhs = [exact(t(1));f(t(2:end-1));exact(t(end))];
for ep=epvec
    A = rbf(ep,r);
    D2 = d2rbf(ep,r);
    Aeval = rbf(ep,reval);
    D2([1,end],:) = A([1,end],:);
    coef = D2\rhs;
    err_Kansa(k) = errcompute(Aeval*coef,uu);
    k = k+1;
end

% Solve the GQRr collocation problem
GQR = gqr_rsolve(t,exact(t),1,1); % Setup GQR
err_GQR = [];
% err_GQR_int = [];
alpha = 1;
rhs = [exact(t(1));f(t(2:end-1));exact(t(end))];
k = 1;
for ep = epvec
%     GQR = gqr_rsolve(t,exact(t),ep,alpha);
%     err_GQR_int(k) = errcompute(gqr_eval(GQR,xx),uu);
    
    [ep,alpha,Marr] = gqr_solveprep(1,t,ep,alpha);
    GQR.ep = ep;
    phiMat = gqr_phi(Marr,t([1,end]),ep,alpha);
    phiMatD2 = gqr_phi(Marr,t(2:end-1),ep,alpha,2);
    
    A = [phiMat(1,:);phiMatD2;phiMat(end,:)];
    GQR.coef = A\rhs;

    err_GQR(k) = errcompute(gqr_eval(GQR,xx),uu);
    k = k+1;
end

warning on MATLAB:nearlySingularMatrix

% Plot the results
loglog(epvec,err_GQR,'r','LineWidth',3),hold on
% loglog(epvec,err_GQR_int,'g-.','LineWidth',2)
loglog(epvec,err_Kansa,'b','LineWidth',2)
loglog(epvec,err_Trefethen*ones(size(epvec)),'--k','LineWidth',2)
hold off
xlabel('\epsilon')
ylabel('absolute 2-norm error')
ylim([1e-13 1e8])
% title(sprintf('Collocation for a 2-pt BVP, N=%d',N))
legend(sprintf('GaussQRr (M=%d)',size(GQR.Marr,2)),'Direct Gauss Collocation','Polynomial Collocation','Location','NorthEast')

figure

rhs = [exact(t(1));f(t(2:end-1));exact(t(end))];
u = D2\rhs;
[pcoef,S,mu] = polyfit(t,u,N-1);
uT = polyval(pcoef,xx,S,mu);

ep = 1;
[ep,alpha,Marr] = gqr_solveprep(1,t,ep,alpha);
GQR.ep = ep;
phiMat = gqr_phi(Marr,t([1,end]),ep,alpha);
phiMatD2 = gqr_phi(Marr,t(2:end-1),ep,alpha,2);
A = [phiMat(1,:);phiMatD2;phiMat(end,:)];
GQR.coef = A\rhs;
uG = gqr_eval(GQR,xx);

% [AX,H1,H2] = plotyy(xx,exact(xx),xx,[abs(uu-uT),abs(uu-uG)],'plot','semilogy');
% set(get(AX(1),'ylabel'),'String','True Solution')
% set(get(AX(2),'ylabel'),'String','Solution Error')
% set(H1,'linewidth',2)
% set(H2,'linewidth',2)
semilogy(xx,abs(uu-uT),'k','linewidth',2),hold on
semilogy(xx,abs(uu-uG)+eps,'r','linewidth',2),hold off
ylabel('pointwise solution error')
ylim([1e-16,1e-6])
xlabel('x')
legend('Polynomial collocation','GaussQRr (M=40,\epsilon=1)','location','north')
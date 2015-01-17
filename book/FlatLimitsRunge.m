% FlatLimitsRunge
% This example demonstrates the ability of ep>0 to alleviate the Runge
% phenomenom even on evenly distributed points
% All the gqr solves use alpha=1
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

% Define the troublesome Runge function
yf = @(x) 1./(1+25*x.^2);

% Define some points on which to interpolate
N = 15;
x = pickpoints(-1,1,N);
y = yf(x);
xcheb = pickpoints(-1,1,N,'cheb');
ycheb = yf(xcheb);

% Define points at which to evaluate and judge the error
Neval = 500;
xeval = pickpoints(-1,1,Neval);
yeval = yf(xeval);

% Plot the polynomial interpolants on different node sets
ypint = polyval(polyfit(x,y,N-1),xeval);
ypint_cheb = polyval(polyfit(xcheb,ycheb,N-1),xeval);
h_poly = figure;
hold on
plot(xcheb,ycheb,'+k','linewidth',1,'markersize',13)
plot(x,y,'ok','linewidth',1,'markersize',13)
plot(xeval,ypint_cheb,'m','linewidth',3)
plot(xeval,ypint,'--','linewidth',3)
hold off
legend('Chebyshev data','uniform data','Chebyshev interpolant',...
       'uniform interpolant','location','north')

% Create Gaussian interpolants with various degrees of locality
% For small ep (near the polynomial limit), this requires GaussQR
ykp1 = gqr_eval(gqr_solve(x,y,.1,1),xeval);
yk1 = gqr_eval(gqr_solve(x,y,1,1),xeval);
yk3 = gqr_eval(gqr_solve(x,y,3,1),xeval);

% Plot the different kernel interpolants together to show that a more
% localized kernel can produce Chebyshev points-like results
h_kernels = figure;
hold on
plot(xeval,yk3,'-g','linewidth',3)
plot(xeval,yk1,'--b','linewidth',3)
plot(xeval,ykp1,':r','linewidth',3)
plot(x,y,'ok','linewidth',1,'markersize',13)
hold off
legend('\epsilon=3 interpolant','\epsilon=1 interpolant',...
       '\epsilon=.1 interpolant','uniform data','location','north')

% Choose some larger epsilon points and solve with the direct method
DM = DistanceMatrix(x,x);
DMeval = DistanceMatrix(xeval,x);
epvecD = logspace(log10(3),2,160);
errvecD = arrayfun(@(ep) errcompute(rbf(ep,DMeval)*(rbf(ep,DM)\y),yeval),epvecD);
   
% Create an epsilon profile, showing the error in the interpolant as a
% function of epsilon
epvecHS = logspace(-1,log10(3),40);
errvecHS = arrayfun(@(ep)errcompute(gqr_eval(gqr_solve(x,y,ep,1),xeval),yeval),epvecHS);

% Compose the direct and stable results and plot them
epvec = [epvecHS,epvecD];
errvec = [errvecHS,errvecD];
h_err = figure;
loglog(epvec,errvec,'linewidth',3)
hold on
loglog(epvec,ones(size(epvec))*errcompute(ypint,yeval),':k')
loglog(epvec,ones(size(epvec))*errcompute(ypint_cheb,yeval),'--k')
hold off
xlabel('\epsilon')
ylabel('absolute 2-norm error')
legend('uniform Gaussian','uniform poly','Chebyshev poly','location','southwest')
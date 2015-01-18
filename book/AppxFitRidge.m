% AppxFitRidge.m
% This program demonstrates Ridge Regression as interpreted two ways
% It is drawn from Grace Wahba's 1990 book, Chapter 4.1

% Set the random number generator to produce consistent results
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.RANDOM_SEED(0);
GAUSSQR_PARAMETERS.ERROR_STYLE = 4;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

% Set up the data in the problem
N = 100;
yf = @(x) 4.26*(exp(-x)-4*exp(-2*x)+3*exp(-3*x));
sigma = .2;
I = eye(N);

x = pickpoints(0,3.5,N);
y = yf(x) + sigma^2*randn(N,1);

% Set up evaluation points at which to study the answer
NN = 400;
xeval = pickpoints(0,3.5,NN);
yeval = yf(xeval);

% Set up kernels with which to approximate
% Also choose a smoothing parameter
ep = .4;
mu = 1e-5;
rbfM2 = @(e,r) (1+e*r).*exp(-(e*r));
rbfM4 = @(e,r) (1+e*r+(e*r).^2/3).*exp(-(e*r));
rbfG = @(e,r) exp(-(e*r).^2);

% Choose a kernel for the approximation
rbf = rbfM4;

% Create the interpolation matrix
DM = DistanceMatrix(x,x);
K = rbf(ep,DM);

% Create the evaluation matrix
DMeval = DistanceMatrix(xeval,x);
Keval = rbf(ep,DMeval);

% Form the standard kernel interpolant
yint = Keval*(K\y);

% Evaluate the standard interpolant and the smoothing spline
% (K+mu*I) is the smoothing spline matrix
Kmu = K+mu*I;
yridge = Keval*(Kmu\y);

% Alternately we can consider the true interpolant to the smoothed data
% We create the smoothed data by noticing that the Smoothing spline is
%     s(x) = u(x)*inv(K+mu*I)*K*y
% So the smoothed data is inv(K+mu*I)*K*y instead of the noisy data y
ysmooth = Kmu\(K*y);

% Plot the data points
h = figure;
subplot(1,3,1)
plot(xeval,yf(xeval),'k','linewidth',3)
hold on
plot(x,ysmooth,'xb','linewidth',3)
plot(x,y,'or','linewidth',1.5)
hold off
xlim([min(x),max(x)])
ylim([-1.2,.4])
legend('true function','smoothed data','noisy data','location','southeast')
title('noisy data')

% Plot the noisy data and the fit
subplot(1,3,2)
plot(x,y,'or')
hold on
plot(xeval,yridge,'linewidth',3)
hold off
xlim([min(x),max(x)])
ylim([-1.2,.4])
title(sprintf('Smoothing spline to noisy data, mu=%g',mu))

% Plot the interpolation to the smoothed data
subplot(1,3,3)
plot(x,ysmooth,'or','markersize',12)
hold on
plot(xeval,yridge,'linewidth',3)
plot(xeval,yeval,'--','linewidth',3)
hold off
xlim([min(x),max(x)])
ylim([-1.2,.4])
legend('smoothed data','ridge regression','true function','location','southeast')
title(sprintf('Interpolant to smoothed data, mu=%g',mu))

% Consider the ridge regression over a variety of ep and mu parameters
epvec = logspace(-3,1,30);
muvec = logspace(-15,0,4);
[E,M] = meshgrid(epvec,muvec);
emvec = num2cell([E(:),M(:)],2);
errvecem = cellfun(@(em)errcompute(rbf(em(1),DMeval)*((rbf(em(1),DM)+em(2)*I)\y),yeval),emvec);
ERR = reshape(errvecem,length(muvec),length(epvec));
h_surf = figure;
surf(epvec,muvec,ERR)
set(gca,'xscale','log','yscale','log','zscale','log')
xlabel('\epsilon'),ylabel('\mu'),zlabel('2-norm error')
set(gca,'xtick',[.01,1],'ytick',[1e-14,1e-7,1])
colormap gray,colormap(flipud(colormap))
view([-.8,-1.2,1])

[epmuopt,besterror] = fminsearch(@(em)errcompute( ...
    rbf(exp(em(1)),DMeval)*((rbf(exp(em(1)),DM)+exp(em(2))*I)\y), ...
    yeval),[-1,-7])
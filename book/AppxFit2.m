% AppxFit2.m
% This program is a basic example which involves Ridge Regression
% It is drawn from Grace Wahba's book, Chapter 4.1

% Set the random number generator to produce consistent results
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.RANDOM_SEED(0);

% Set up the data in the problem
N = 100;
yf = @(x) 4.26*(exp(-x)-4*exp(-2*x)+3*exp(-3*x));
sigma = .2;
I = eye(N);

x = pickpoints(0,3.5,N);
y = yf(x) + sigma^2*randn(N,1);

% Set up evaluation points at which to study the answer
NN = 400;
xx = pickpoints(0,3.5,NN);
yy = yf(xx);

% Set up kernels with which to approximate
% Also choose a smoothing parameter
ep = 1;
mu = 1e-3;
rbfM2 = @(e,r) (1+e*r).*exp(-(e*r));
rbfG = @(e,r) exp(-(e*r).^2);

% Choose a kernel for the approximation
rbf = rbfG;

% Create the interpolation matrix
DM = DistanceMatrix(x,x);
K = rbf(ep,DM);

% Create the evaluation matrix
DM_eval = DistanceMatrix(xx,x);
K_eval = rbf(ep,DM_eval);

% Evaluate the smoothing spline
% (K+mu*I) is the smoothing spline matrix
K_mu = K+mu*I;
y_eval = K_eval*(K_mu\y);

% Alternately we can consider the true interpolant to the smoothed data
% We create the smoothed data by noticing that the Smoothing spline is
%     s(x) = u(x)*inv(K+mu*I)*K*y
% So the smoothed data is inv(K+mu*I)*K*y instead of the noisy data y
y_smooth = K_mu\(K*y);

% If we use the smoothed out data, we can consider implementing the stable
% basis to evaluate the cardinal functions which would allow us to evaluate
% the interpolant
% Here, we are going to evaluate the cardinal functions as though the data
% were shifted, but not scaled, to be centered around 0
% cardinal_functions contains the cardinal functions at the evaluation
% points
% x_zeroed = (x - 3.5/2)/(3.5/2);
% xx_zeroed = (xx - 3.5/2)/(3.5/2);
% gqr_alpha = 1.01;
% gqr_M = 25;
% GQR = gqr_solveprep(0,x_zeroed,ep,gqr_alpha);
% Psi = GQR.stored_phi1 + GQR.stored_phi2*GQR.Rbar;
% Phi_eval = gqr_phi(GQR,x_zeroed);
% Psi_eval = Phi_eval*[I;GQR.Rbar];
% cardinal_functions = Psi_eval/Psi;
% GQR = gqr_solveprep(1,x_zeroed,ep,gqr_alpha,gqr_M);
% Phi = gqr_phi(GQR,x_zeroed);
% Phi_eval = gqr_phi(GQR,xx_zeroed);
% cardinal_functions = Phi_eval/Phi;
% figure,plot(xx,cardinal_functions(:,[1,6,28,62,99])),hold on,plot(x,zeros(size(x)),'or'),hold off

% Plot the data points
h = figure;
subplot(1,3,1)
plot(xx,yf(xx),'k','linewidth',3)
hold on
plot(x,y,'or','linewidth',1.5)
plot(x,y_smooth,'xb','linewidth',3)
hold off
title('noisy data')
xlim([min(x),max(x)])
ylim([-1.2,.4])

% Plot the noisy data and the fit
subplot(1,3,2)
plot(x,y,'or')
hold on
plot(xx,y_eval,'linewidth',3)
hold off
title(sprintf('Smoothing spline to noisy data, mu=%g',mu))
xlim([min(x),max(x)])
ylim([-1.2,.4])

% Plot the interpolation to the smoothed data
subplot(1,3,3)
plot(x,y_smooth,'or')
hold on
plot(xx,y_eval,'linewidth',3)
hold off
title(sprintf('Interpolant to smoothed data, mu=%g',mu))
xlim([min(x),max(x)])
ylim([-1.2,.4])
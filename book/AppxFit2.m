% AppxFit2.m
% This program is a basic example which involves Ridge Regression
% It is drawn from Grace Wahba's book, Chapter 4.1
N = 100;
yf = @(x) 4.26*(exp(-x)-4*exp(-2*x)+3*exp(-3*x));
sigma = .2;
I = eye(N);

x = pickpoints(0,3.5,N);
y = yf(x) + sigma^2*randn(N,1);

NN = 400;
xx = pickpoints(0,3.5,NN);
yy = yf(xx);

% Choose a kernel with which to approximate
% Also choose a smoothing parameter
ep = 1;
mu = 1e-1;
rbfM2 = @(e,r) (1+e*r).*exp(-(e*r));

% Create the interpolation matrix
DM = DistanceMatrix(x,x);
K = rbfM2(ep,DM);

% Create the evaluation matrix
DM_eval = DistanceMatrix(xx,x);
K_eval = rbfM2(ep,DM_eval);

% Evaluate the smoothing spline
% (K+mu*I) is the smoothing spline matrix
K_mu = K+mu*I;
y_eval = K_eval*(K_mu\y);

% Alternately we can consider the true interpolant to the smoothed data
% We create the smoothed data by noticing that the Smoothing spline is
%     s(x) = u(x)*inv(K+mu*I)*K*y
% So the smoothed data is inv(K+mu*I)*K*y instead of the noisy data y
y_smooth = K_mu\(K*y);

% Plot the data points
h = figure;
subplot(1,3,1)
plot(xx,f(xx),'k','linewidth',3)
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
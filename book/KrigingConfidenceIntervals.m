% KrigingConfidenceIntervals
% This example demonstrates the use of confidence intervals to show how
% much you trust your Kriging predictions

% Choose a kernel, here the Gaussian
rbf = @(e,r) exp(-(e*r).^2);
ep = 3;
sigma = 1;

% Points at which to sample
N = 10;
x = pickpoints(-1,1,N,'cheb');
Neval = 100;
xeval = unique([pickpoints(-1,1,Neval);x]);

% Choose a function and sample it
yf = @(x) exp(x) - cos(2*pi*x);
y = yf(x);

% Compute the interpolant through those points
DM = DistanceMatrix(x,x);
DMeval = DistanceMatrix(xeval,x);
K = rbf(ep,DM);
Keval = rbf(ep,DMeval);
seval = Keval*(K\y);

% Compute the standard deviation at all the evaluation locations
sd = sigma*real(sqrt(rbf(ep,0) - sum((Keval/K).*Keval,2)));

% Plot the results
h = figure;
plot(x,y,'or');
hold on
plot(xeval,seval,'linewidth',1);
plot(xeval,seval+2*sd,'--k','linewidth',2);
plot(xeval,seval-2*sd,'--k','linewidth',2);
hold off
legend('Data','Prediction','+/-2 SD','location','southeast')
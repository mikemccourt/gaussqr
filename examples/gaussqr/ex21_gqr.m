% ex21_gqr.m
% This example demonstrates the use of confidence intervals to show how
% much you trust your Kriging predictions

% Choose an RBF
% I might wanna try coding this up to not use RBFs only
rbf = @(e,r) exp(-(e*r).^2);
ep = 3;

% Points at which to sample
N = 10;
x = pickpoints(-1,1,N,'cheb');
NN = 100;
xx = unique([pickpoints(-1,1,NN);x]);

% Choose a function and sample it
yf = @(x) exp(x) - cos(2*pi*x);
y = yf(x);
yy = yf(xx);

% Compute the interpolant through those points
DM_INT = DistanceMatrix(x,x);
DM_EVAL = DistanceMatrix(xx,x);
IM = rbf(ep,DM_INT);
EM = rbf(ep,DM_EVAL);
yp = EM*(IM\y);

% Compute the standard deviation at all the evaluation locations
sd = sqrt(rbf(ep,0) - sum((EM/IM).*EM,2));

% Plot the results
h = figure;
plot(x,y,'or');
hold on
plot(xx,yp,'linewidth',1);
plot(xx,yp+2*sd,'--k','linewidth',2);
plot(xx,yp-2*sd,'--k','linewidth',2);
hold off
legend('Data','Prediction','+/-2 SD','location','southeast')
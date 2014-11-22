% StatFit3.m
% This uses the carsmall dataset in Matlab to conduct some statistical data
% fitting on the impact of:
%     [Acceleration Displacement Horsepower Weight]
% on a car's MPG.

% Load, clean and scale the data
load carsmall
xdirty = [Acceleration Displacement Horsepower Weight];
ydirty = MPG;
[x,y,shift,scale] = rescale_data(xdirty,ydirty);

% Try the interpolation in just 2D for now
NN = 50;
xx = pick2Dpoints([-1 -1],[1 1],NN*ones(1,2));
ep = 5;
rbf = @(e,r) exp(-(e*r).^2);
[x_AD,ind_AD] = unique(x(:,1:2),'rows');

N_AD = length(ind_AD);
y_AD = y(ind_AD);
K = rbf(ep,DistanceMatrix(x_AD,x_AD));
mu = 1e-3;

yp = rbf(ep,DistanceMatrix(xx,x_AD))*((K'*K+mu*eye(N_AD))\(K'*y_AD));
h_AD = figure;
surf(reshape(xx(:,1),NN,NN),reshape(xx(:,2),NN,NN),reshape(yp,NN,NN))


alpha = 1;
GQR_AD = gqr_solve(x_AD,y_AD,ep,alpha);
yp = gqr_eval(GQR_AD,xx);
h_AD = figure;
surf(reshape(xx(:,1),NN,NN),reshape(xx(:,2),NN,NN),reshape(yp,NN,NN))
% SurrModelCarSmall
% This example considers data from Matlab's carsmall data set containing
% the MPG ratings for cars with various parameters:
%       {'Acceleration','Displacement','Horsepower','Weight'}
% We conduct some surrogate modeling and study plots of the model

% Define some radial kernels for use on this problem
rbfM2 = @(r) (1+r).*exp(-r);
rbfM2dx = @(r,dx,ep) -ep^2*exp(-r).*dx;
rbfM2dxdy = @(r,dx,dy,ep) prod(ep.^2)*exp(-r).*dx.*dy./(r+eps);
rbfM4 = @(r) (3+3*r+r.^2).*exp(-r);
rbfM4dx = @(r,dx,ep) -ep^2*exp(-r).*(1+r).*dx;
rbfM4dxdy = @(r,dx,dy,ep) prod(ep.^2)*exp(-r).*dx.*dy;
rbfM6 = @(r) (15+15*r+6*r.^2+r.^3).*exp(-r);
rbfM6dx = @(r,dx,ep) -ep^2*exp(-r).*(r.^2+3*r+3).*dx;
rbfM6dxdy = @(r,dx,dy,ep) prod(ep.^2)*exp(-r).*dx.*dy.*(1+r);
rbfG = @(r) exp(-r.^2);

% This function allows you to evaluate the EDF
% Here, xe are the evaluation points, x are the observed locations
Fhat = @(xe,x) reshape(sum(all(repmat(x,[1,1,size(xe,1)])<=repmat(reshape(xe',[1,size(xe,2),size(xe,1)]),[size(x,1),1,1]),2),1),size(xe,1),1)/size(x,1);

% Load, clean and scale the data
load carsmall
xdirty = [Acceleration Displacement Horsepower Weight];
xstr = {'Acceleration','Displacement','Horsepower','Weight'};
ydirty = MPG;
[x,y,shift,scale] = rescale_data(xdirty,ydirty);
x_mean = mean(x);

% Try to compute a 2D density over acceleration and horsepower
xAD = x(:,1:2);
NAD = size(xAD,1);
N2d = 35;
x2d = pick2Dpoints([-1 -1],[1 1],N2d*[1;1]);

% Sorting may be useful but I haven't figured out why yet
[~,i] = sort(sum(x2d - ones(N2d^2,1)*[-1,-1],2));
x2d_sorted = x2d(i,:);
[c,i] = sort(sum(xAD - ones(NAD,1)*[-1,-1],2));
xAD_sorted = xAD(i,:);

% Create a scatter plot of the data locations in 2D
h_scatter = figure;
scatter(xAD_sorted(:,1),xAD_sorted(:,2),exp(3*c))
ecdf2d = zeros(N2d^2,1);
for k=1:N2d^2
    ecdf2d(k) = sum(all(xAD<=repmat(x2d(k,:),NAD,1),2))/NAD;
end
h_ecdf = figure;
surf(reshape(x2d(:,1),N2d,N2d),reshape(x2d(:,2),N2d,N2d),reshape(ecdf2d,N2d,N2d))

% Choose a radial kernel
rbf = rbfM6;
rbfdxdy = rbfM6dxdy;

% Compute a smoothing spline surrogate model on the distribution of the
% data in 2D, and evaluate the PDF
ep = [3,3];
mu = 1e-3;
K_cdf = rbf(DistanceMatrix(x2d,x2d,ep));
cdf2d_coef = (K_cdf+mu*eye(N2d^2))\ecdf2d;
cdf2d_eval = @(xx) rbf(DistanceMatrix(xx,x2d,ep))*cdf2d_coef;
pdf2d_eval = @(xx) max(rbfdxdy(DistanceMatrix(xx,x2d,ep),DifferenceMatrix(xx(:,1),x2d(:,1)),...
    DifferenceMatrix(xx(:,2),x2d(:,2)),ep)*cdf2d_coef,0);
Neval = 50;
x2d_eval = pick2Dpoints([-1 -1],[1 1],Neval*[1;1]);
y_eval = cdf2d_eval(x2d_eval);
surf(reshape(x2d_eval(:,1),Neval,Neval),reshape(x2d_eval(:,2),Neval,Neval),reshape(y_eval,Neval,Neval))
y_eval = pdf2d_eval(x2d_eval);
surf(reshape(x2d_eval(:,1),Neval,Neval),reshape(x2d_eval(:,2),Neval,Neval),reshape(y_eval,Neval,Neval))

%%%%%%%%%%
% Shift to using the Gaussians to consider the relationship between
% dimensions in the full 4D surrogate model
rbf = rbfG;

% Create the surrogate model using a pre-chosen shape parameter
ep = 4;
K = rbf(DistanceMatrix(x,x,ep));
coef = K\y;
SM_eval = @(xx) rbf(DistanceMatrix(xx,x,ep))*coef;

% Evaluate cross-sections of the surrogate model
% Two of the inputs will be fixed to their means
% The other two will be allowed to vary over their input range
NN = 50;
x2d = pick2Dpoints([-1 -1],[1 1],NN*[1;1]);
X1 = reshape(x2d(:,1),NN,NN);
X2 = reshape(x2d(:,2),NN,NN);
mean_vecs = ones(NN^2,1)*x_mean;
subplot_cs = {[1 2],[1 3],[1 4],[2 3],[2 4],[3 4]};
h_cs = figure;
for k=1:length(subplot_cs)
    subplot(3,2,k)
    this_subplot = subplot_cs{k};
    xeval = mean_vecs;
    xeval(:,this_subplot) = x2d;
    yeval = SM_eval(xeval);
    surf(X1,X2,reshape(yeval,NN,NN),'edgecolor','none');
    xlabel(xstr{this_subplot(1)});
    ylabel(xstr{this_subplot(2)});
    zlabel('MPG');
    zlim([0 30])
end

% Define labels for the parallel coordinate plots
xstr = {'Acc','Disp','HP','Weight'};
% Create some parallel coordinate plots
N = 15;
x1 = pickpoints(-1,1,N);
[X1,X2,X3,X4] = ndgrid(x1);
x4d = [X1(:),X2(:),X3(:),X4(:)];
yfull = SM_eval(x4d);
subplot_pc = {[9 14],[14 21],[21 28],[28 55]};
h_pc = figure;
for k=1:length(subplot_pc)
    subplot(2,2,k);
    this_subplot = subplot_pc{k};
    this_ind = yfull>this_subplot(1) & yfull<this_subplot(2);
    xgood = x4d(this_ind,:);
    groups = ceil(3*(max(xgood(:,1)+1,eps))./2);
    parallelcoords(xgood,'Group',groups,'Label',xstr);
    title(sprintf('MPG between %d and %d',this_subplot))
end
% SurrModel4.m
% This considers the development of a surrogate model for the empirical
% distribution function (ECDF) of data drawn from a random distribution.
% The data we consider again comes from the carsmall data set, where we
% consider the probability of finding a car with a given set of
% Acceleration, Displacement, Horsepower and Weight values.
% We use this density to help us compute a marginal model, studying the
% relevance of the parameters Horsepower and Weight averaged over all
% possible Acceleration and Displacement values which occur with
% probability density defined through the ECDF surrogate model.

N = 200;
x = randn(N,2);
xe = pick2Dpoints(min(x(:))*[1 1],max(x(:))*[1 1],[30;30]);
Fhat = @(xe,x) reshape(sum(all(repmat(x,[1,1,size(xe,1)])<=repmat(reshape(xe',[1,2,size(xe,1)]),[size(x,1),1,1]),2),1),size(xe,1),1)/size(x,1);

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
[c,i] = sort(sum(x2d - ones(N2d^2,1)*[-1,-1],2));
x2d_sorted = x2d(i,:);
[c,i] = sort(sum(xAD - ones(NAD,1)*[-1,-1],2));
xAD_sorted = xAD(i,:);
h_scatter = figure;
scatter(xAD_sorted(:,1),xAD_sorted(:,2),exp(3*c))
ecdf2d = zeros(N2d^2,1);
for k=1:N2d^2
    ecdf2d(k) = sum(all(xAD<=repmat(x2d(k,:),NAD,1),2))/NAD;
end
h_ecdf = figure;
surf(reshape(x2d(:,1),N2d,N2d),reshape(x2d(:,2),N2d,N2d),reshape(ecdf2d,N2d,N2d))

rbfM2 = @(r) (1+r).*exp(-r);
rbfM2dxdy = @(r,dx,dy,ep) prod(ep.^2)*exp(-r).*dx.*dy./(r+eps);
rbfM4 = @(r) (3+3*r+r.^2).*exp(-r);
rbfM4dxdy = @(r,dx,dy,ep) prod(ep.^2)*exp(-r).*dx.*dy;
rbfM6 = @(r) (15+15*r+6*r.^2+r.^3).*exp(-r);
rbfM6dxdy = @(r,dx,dy,ep) prod(ep.^2)*exp(-r).*dx.*dy.*(1+r);
rbf = rbfM6;
rbfdxdy = rbfM6dxdy;

ep = [3,3];
K_cdf = rbf(DistanceMatrix(x2d,x2d,ep));
cdf2d_coef = K_cdf\ecdf2d;
cdf2d_eval = @(xx) rbf(DistanceMatrix(xx,x2d,ep))*cdf2d_coef;
pdf2d_eval = @(xx) max(rbfdxdy(DistanceMatrix(xx,x2d,ep),DifferenceMatrix(xx(:,1),x2d(:,1)),...
                             DifferenceMatrix(xx(:,2),x2d(:,2)),ep)*cdf2d_coef,0);
Neval = 50;
x2d_eval = pick2Dpoints([-1 -1],[1 1],Neval*[1;1]);
y_eval = cdf2d_eval(x2d_eval);
surf(reshape(x2d_eval(:,1),Neval,Neval),reshape(x2d_eval(:,2),Neval,Neval),reshape(y_eval,Neval,Neval))
y_eval = pdf2d_eval(x2d_eval);
surf(reshape(x2d_eval(:,1),Neval,Neval),reshape(x2d_eval(:,2),Neval,Neval),reshape(y_eval,Neval,Neval))
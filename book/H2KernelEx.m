f = @(x) cos(x)+exp(-(x-1).^2)-exp(-(x+1).^2);
N = 10;
x = pickpoints(-2,2,N);
y = f(x);
NN = 150;
xx = pickpoints(-2,2,NN);
yy = f(xx);

rbf1 = @(ep,r) exp(-ep*r).*(1+ep*r);    % C^2 Matern
rbf2 = @(ep,r) sqrt(2)*exp(-ep*r).*sin(ep*r + pi/4);    % H^2 Sobolev spline (BTA)
ep = 3;

DM = DistanceMatrix(x,x);
DMeval = DistanceMatrix(xx,x);

figure
K = rbf1(ep,DM);
Keval = rbf1(ep,DMeval);
yp = Keval*(K\y);
hold on
plot(xx,yy,'linewidth',3)
plot(xx,yp,'linewidth',3)
hold off
title(sprintf('Error Matern = %g',errcompute(yp,yy)))

figure
K = rbf2(ep,DM);
Keval = rbf2(ep,DMeval);
yp = Keval*(K\y);
hold on
plot(xx,yy,'linewidth',3)
plot(xx,yp,'linewidth',3)
hold off
title(sprintf('Error Sobolev = %g',errcompute(yp,yy)))

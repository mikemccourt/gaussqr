% FlatLimitsRunge
% This example demonstrates the ability of ep>0 to alleviate the Runge
% phenomenom

yf = @(x) 1./(1+25*x.^2);

N = 15;
x = pickpoints(-1,1,N);
y = yf(x);

Neval = 500;
xx = pickpoints(-1,1,Neval);
yy = yf(xx);

ypp1 = gqr_eval(gqr_solve(x,y,.1,1),xx);
yp1 = gqr_eval(gqr_solve(x,y,1,1),xx);
yp3 = gqr_eval(gqr_solve(x,y,3,1),xx);

h = figure;
hold on
plot(xx,yy,'linewidth',3)
plot(xx,[ypp1,yp1],'linewidth',3)
plot(xx,yp3,'--m','linewidth',3)
plot(x,y,'ok','linewidth',1,'markersize',13)
hold off
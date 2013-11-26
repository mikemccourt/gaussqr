%Computes the eigenvalues of K(x,z) = exp(-a*|x-z|) on [-1,1] by solving
%transcendental equations from Dongbin Xiu, Chapter 4 (eqns. 4.9 and 4.10)
clear all
close all
nmax = 20;
w = zeros(1,nmax);
a = 1;

% Even indices
ycrop = a*nmax*pi;
fun = @(x) tan(x)+a*x;
xx = linspace(0,(2*nmax+1)*pi/2,500);
hold on
yy = tan(xx);
jump = find(abs(yy)>ycrop);
yy(jump) = NaN;
plot(xx,yy,'r','LineWidth',2)
text(0.9*(2*nmax+1)*pi/2,-.5,'w','fontsize',14)
set(gca,'XGrid','on')
set(gca,'YGrid','on')
axis([0  (2*nmax+1)*pi/2  -ycrop  ycrop])
plot(xx,-a*xx,'r','LineWidth',2)
% Let's compute the eigenvalues numerically
disp(sprintf(' n       w       (2n-1)pi/2    lambda(n)'))
for n=1:nmax
    w(n) = fzero(fun,[(2*n-1)*pi/2+100*eps (2*n+1)*pi/2-100*eps]);
    disp(sprintf('%2d  %10f   %10f  %10f',n, w(n), (2*n-1)*pi/2, 2*a/(1+a*a*w(n)*w(n))))
end
plot(w,-a*w,'bo','LineWidth',2)
hold off

pause
figure
v = zeros(1,nmax+1);
% Odd indices
ycrop = 5;
fun = @(x) 1-a*x*tan(x);
xx = linspace(0,(2*nmax+1)*pi/2,500);
hold on
yy = tan(xx);
jump = find(abs(yy)>ycrop);
yy(jump) = NaN;
plot(xx,yy,'r','LineWidth',2)
text(0.95*(2*nmax+1)*pi/2,-.5,'v','fontsize',14)
set(gca,'XGrid','on')
set(gca,'YGrid','on')
axis([0  (2*nmax+1)*pi/2  -1  ycrop])
plot(xx,1./(a*xx),'r','LineWidth',2)
% Let's compute the eigenvalues numerically
disp(sprintf(' n       v           n*pi      lambda(n)'))
v(1) = fzero(fun,[0 *pi/2-100*eps]);
disp(sprintf('%2d  %10f   %10f  %10f',1, v(1), 0, 2*a/(1+a*a*v(1)*v(1))))
for n=1:nmax
    v(n+1) = fzero(fun,[(2*n-1)*pi/2+100*eps (2*n+1)*pi/2-100*eps]);
    disp(sprintf('%2d  %10f   %10f  %10f',n+1, v(n+1), n*pi, 2*a/(1+a*a*v(n+1)*v(n+1))))
end
plot(v,1./(a*v),'bo','LineWidth',2)
hold off

figure
lambda = zeros(1,2*nmax+1);
for n=1:nmax
    lambda(2*n-1) = 2*a/(1+a*a*v(n)*v(n));
    lambda(2*n) = 2*a/(1+a*a*w(n)*w(n));
end
lambda(2*nmax+1) = 2*a/(1+a*a*v(nmax+1)*v(nmax+1));
semilogy(1:(2*nmax+1),lambda)
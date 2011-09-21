% ex7
% This example produces the surface graph from Chapter 6.2

yf = @(x) 10*exp(-x.^2)+x.^2;
N = 200;
x = pickpoints(-5,5,N);
y = yf(x);
NN = 1000;
xx = pickpoints(-5,5,NN);
yy = yf(xx);

alpha = 1;
epvec = logspace(-1,1,40);
Mvec = round(linspace(20,180,39));
errs = zeros(length(epvec),length(Mvec));

i = 1;
for ep=epvec
    j = 1;
    for M=Mvec
        RBF = rbfqrr_solve_alpha(x,y,ep,alpha,M);
        yp = rbfqr_eval_alpha(RBF,xx);
        errs(i,j) = norm((yy-yp)./(abs(yp)+eps))/NN;
        j = j+1;
    end
    fprintf('%d ',i)
    i = i+1;
end

figure
[EE,MM] = meshgrid(epvec,Mvec);
surfc(EE,MM,log10(errs)')
xlabel('\epsilon')
ylabel('M')
zlabel('log_{10}(error)')
set(gca,'XScale','log')
title('f(x)=10e^{-x^2}+x^2, N=200, \alpha=1, x=[-5,5]')

figure
errC = errs;
errC(find(errs>1)) = ones(size(find(errs>1)));
contour(EE,MM,log10(errC)')
xlabel('\epsilon')
ylabel('M')
set(gca,'XScale','log')
title('f(x)=10e^{-x^2}+x^2, N=200, \alpha=1, x=[-5,5]')
h = colorbar;
set(get(h,'Ylabel'),'String','log_{10}(error)')
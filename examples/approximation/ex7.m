% ex7
% This example studies the relationship between epsilon and the
% eigenfunction basis on approximation quality

yf = @(x) 10*exp(-x.^2)+x.^2;
N = 200;
x = pickpoints(-5,5,N);
y = yf(x);
Neval = 1000;
xeval = pickpoints(-5,5,Neval);
yeval = yf(xeval);

alpha = 1;
epvec = logspace(-1,1,20);
Mvec = round(linspace(20,180,19));
errs = zeros(length(epvec),length(Mvec));

h_waitbar = waitbar(0,'Initializing');
i = 1;
for ep=epvec
    j = 1;
    for M=Mvec
        RBF = gqr_rsolve(x,y,ep,alpha,M);
        ypred = gqr_eval(RBF,xeval);
        errs(i,j) = errcompute(ypred,yeval);
        waitbar((j+(i-1)*length(Mvec))/(length(Mvec)*length(epvec)),h_waitbar,sprintf('ep=%3.2f, M=%d',ep,M))
        j = j+1;
    end
    i = i+1;
end
waitbar(100,h_waitbar,'Plotting')

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
errC(errs>1) = 1;
contour(EE,MM,log10(errC)')
xlabel('\epsilon')
ylabel('M')
set(gca,'XScale','log')
title('f(x)=10e^{-x^2}+x^2, N=200, \alpha=1, x=[-5,5]')
h = colorbar;
set(get(h,'Ylabel'),'String','log_{10}(error)')

close(h_waitbar)
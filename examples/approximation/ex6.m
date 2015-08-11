% ex6
% This is an example showing that regression should be a viable option
% because the underlying function only has some maximum amount of
% complexity that you might be able to capture with limited basis
% functions.
%
% Here a function is considered with several N values of input points, and
% several M values of regression complexity
global GAUSSQR_PARAMETERS
Mextramax = GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC;

h_waitbar = waitbar(0,'Initializing');

Nrange = 10:10:500;
Neval = 200;
ep = 1;
alpha = 1;
yf = @(x) cos(x)+exp(-(x-1).^2)+exp(-(x+1).^2);
Mfrac = .1:.1:.5;
errs = zeros(length(Mfrac),length(Nrange));

for i=1:length(Mfrac)
    for j=1:length(Nrange)
        N = Nrange(j);
        M = round(N*Mfrac(i));
        waitbar(((i-1)*length(Nrange)+j)/(length(Nrange)*length(Mfrac)),h_waitbar,sprintf('N=%d, M=%d',N,M));
        
        x = pickpoints(-3,3,N);
        y = yf(x);
        xeval = pickpoints(-3,3,Neval,'halt');
        yeval = yf(xeval);
        RBFQR = gqr_rsolve(x,y,ep,alpha,M);
        yp = gqr_eval(RBFQR,xeval);
        errs(i,j) = errcompute(yp,yeval);
    end
end

waitbar(1,h_waitbar,'Plotting')

semilogy(Nrange,errs)
title(sprintf('\\alpha=%g, \\epsilon=%g',alpha,ep))
legend('M=.1N','M=.2N','M=.3N','M=.4N','M=.5N','Location','NorthEast')
xlabel('N')
ylabel('Average error')
ylim([1e-16 1])

close(h_waitbar)
% ex6
% This is an example showing that regression should be a viable option
% because the underlying function only has some maximum amount of
% complexity that you might be able to capture with limited basis
% functions.
%
% Here a function is considered with several N values of input points, and
% several M values of regression complexity
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
Mextramax = GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC;

Nrange = 10:10:500;
if not(exist('ep'))
    ep = 1e-1;
end
if not(exist('alpha'))
    alpha = 1;
end
yf = @(x) cos(x)+exp(-(x-1).^2)+exp(-(x+1).^2);
Mfrac = .1:.1:.5;
errs = zeros(length(Mfrac),length(Nrange));

for i=1:length(Mfrac)
    for j=1:length(Nrange)
        N = Nrange(j);
        M = round(N*Mfrac(i));
        x = pickpoints(-3,3,N);
        y = yf(x);
        xx = pickpoints(-3,3,1000);
        yy = yf(xx);
        RBFQR = rbfqrr_solve_alpha(x,y,ep,alpha,M);
        yp = rbfqr_eval_alpha(RBFQR,xx);
        errs(i,j) = norm((yy-yp)./(abs(yp)+eps));
    end
end

semilogy(Nrange,errs)
title(sprintf('\\alpha=%g, \\epsilon=%g',alpha,ep))
legend('M=.1N','M=.2N','M=.3N','M=.4N','M=.5N','Location','NorthEast')
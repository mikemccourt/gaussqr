% ex8
% A study of the condition number of interpolation matrices as a function
% of epsilon and alpha
%
% This example runs slowly, and should be considered only as an outline for
% how to consider ranges of alpha values and that an optimal alpha value
% may exist for conditioning (but not necessarily accuracy purposes)

Nvec = [50 100 200 400];
epvec = [.01 .1 1 10];
alphavec = logspace(-2,2,30);
Mvec = 5:40;
errs = zeros(length(alphavec),length(Mvec));
[AA,MM] = meshgrid(alphavec,Mvec);

k = 1;
for N=Nvec
    x = pickpoints(-5,5,N);
    for ep=epvec
        m = 1;
        for M=Mvec
            a = 1;
            for alpha=alphavec
                Marr = gqr_formMarr(M);
                phi = gqr_phi(Marr,x,ep,alpha);
                errs(a,m) = cond(phi);
                a = a+1;
            end
            m = m+1;
        end
        subplot(length(Nvec),length(epvec),k)
        contour(MM',AA',log10(errs),[2 5 8 11 14],'LineWidth',2)
        xlabel('M')
        ylabel('\alpha')
        set(gca,'Yscale','log')
        title(sprintf('N=%d, \\epsilon=%g',N,ep))
        k = k+1;
    end
end
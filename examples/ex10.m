% Tests orthogonality for different alpha

rbfsetup
weightfunc = @(a,x) a/sqrt(pi)*exp(-a^2*x.^2);
epvec = [.1 1 10]; % Shape parameter
indexvec = [1 5 10]; % Eigenfunction orthogonality to test
Na = 100; % Number of alphas to look at
a = -4;
b = 4;

k=1;
alphavec = logspace(-2,1,Na);
for index=indexvec
    for ep=epvec
        erra = zeros(1,Na);
        j = 1;
        for alpha=alphavec
            erra(j) = quadl(@(x)rbfphialpha(index,x',ep,alpha)'.^2.*weightfunc(alpha,x),a,b);
            j = j+1;
        end
        subplot(length(indexvec),length(epvec),k)
        semilogx(alphavec,erra)
        xlabel('\alpha')
        ylabel('integral appx')
        title(sprintf('ep=%g, index=%d',ep,index))
        k = k+1;
    end
end
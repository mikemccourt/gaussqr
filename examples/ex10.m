% Tests orthogonality for different alpha

rbfsetup
weightfunc = @(a,x) a/sqrt(pi)*exp(-a^2*x.^2);
epvec = [.1 1 10]; % Shape parameter
ptopt = 'even'; % Point distribution to consider
index = 10; % Eigenfunction orthogonality to test
Na = 100; % Number of alphas to look at
a = -3;
b = 3;

k=1;
alphavec = logspace(-2,2,Na);
for ep=epvec
    erra = zeros(1,Na);
    j = 1;
    for alpha=alphavec
        erra(j) = quadl(@(x)rbfphialpha(index,x',ep,alpha)'.*weightfunc(alpha,x).*rbfphialpha(index,x',ep,alpha)',a,b);
        j = j+1;
    end
    subplot(length(epvec),1,k)
    semilogx(alphavec,erra)
    xlabel('\alpha')
    ylabel('integral appx')
    title(sprintf('ep=%g, index=%d',ep,index))
    k = k+1;
end
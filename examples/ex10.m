% Tests orthogonality for different alpha

rbfsetup
weightfunc = @(a,x) a/sqrt(pi)*exp(-a^2*x.^2);
epvec = [.1 1 10]; % Shape parameter
Nvec = [30 100 300]; % Number of points for the trapezoid rule
ptopt = 'even'; % Point distribution to consider
index = 1; % Eigenfunction orthogonality to test
Na = 100; % Number of alphas to look at

k=1;
for N=Nvec
    xx = pickpoints(-3,3,N,ptopt);
    alphavec = logspace(-2,2,Na);
    for ep=epvec
        erra = zeros(1,Na);
        j = 1;
        for alpha=alphavec
            phi_j = rbfphialpha(index,xx,ep,alpha);
            erra(j)=trapz(xx,phi_j.*weightfunc(alpha,xx).*phi_j);
            j = j+1;
        end
        subplot(length(epvec),length(Nvec),k)
        semilogx(alphavec,erra)
        xlabel('\alpha')
        ylabel('integral appx')
        title(sprintf('N=%d, ep=%g, index=%d',N,ep,index))
        k = k+1;
    end
end
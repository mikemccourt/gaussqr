% HSSVDFlatLimit
% This program uses the IBB kernel and tries to compute in the limit as
% ep->0 for a large-ish beta.  The standard basis for the kernel appears
% ill-conditioned so the HS-SVD basis is used for stability
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

% Pick a function that satisfies the IBB boundary conditions
yf = @(x) .25^(-28)*max(x-.25,0).^14.*max(.75-x,0).^14;

% Choose parameters to study
beta = 7;
epvec = logspace(0,3,15);

% Consider different point totals
Nvec = [15,30,60];

% Choose points at which to test the error
Neval = 100;
xeval = pickpoints(0,1,Neval);
yeval = yf(xeval);

% Create the kernel content (for the direct method only)
phifunc = @(n,x) sqrt(2)*sin(pi*x*n);
lamfunc = @(b,e,n) ((pi*n).^2+e^2).^(-b);
M = ceil(1/pi*sqrt(eps^(-1/beta)*(N^2*pi^2+ep^2)-ep^2));

% Compute the error on the N values under consideration
k = 1;
errmathssvd = zeros(length(Nvec),length(epvec));
errmatstandard = zeros(length(Nvec),length(epvec));
for N=Nvec
    x = pickpoints(0,1,N+2);
    x = x(2:end-1);
    y = yf(x);
    
    j = 1;
    for ep=epvec
        Phi = phifunc(1:M,x);
        Phieval = phifunc(1:M,xeval);
        lamvec = lamfunc(beta,ep,1:M);
        K = bsxfun(@times,lamvec,Phi)*Phi';
        Keval = bsxfun(@times,lamvec,Phieval)*Phi';
        errmatstandard(k,j) = errcompute(Keval*(K\y),yeval);
        j = j + 1;
    end
    
    errmathssvd(k,:) = arrayfun(@(ep)errcompute(HSSVD_IBBSolve(ep,beta,x,y,xeval),yeval),epvec);
    
    k = k + 1;
end

% Plot the results on the same plot
h = figure;
h_standard = loglog(epvec,errmatstandard,'--','linewidth',2);
hold on
h_hssvd = loglog(epvec,errmathssvd,'linewidth',2);
hold off
xlabel('$\varepsilon$ - shape parameter','interpreter','latex')
ylabel('2-norm error')
legend([h_standard(1),h_hssvd(1)],{'standard basis','HS-SVD basis'},...
       'location','north')
set([h_hssvd(2),h_standard(2)],'marker','o')
set([h_hssvd(3),h_standard(3)],'marker','+')
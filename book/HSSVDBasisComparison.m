% HSSVDBasisComparison
% This example studies the effectiveness of alternate bases on computing
% interpolants in the flat limit
% We will compare the standard basis, Newton basis, SVD basis,
% eigenfunction basis and HS-SVD basis
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

% Define the Gaussian kernel
rbf = @(e,r) exp(-(e*r).^2);

% Define the Gauss-Hermite quadrature locations and weights
% These are recommended for the weighted SVD basis
GH_int = ...
[-5.3874808900112	2.229393645534E-13
-4.6036824495507	4.399340992273E-10
-3.9447640401156	1.086069370769E-7
-3.3478545673832	7.80255647853E-6
-2.7888060584281	2.283386360164E-4
-2.2549740020893	0.003243773342238
-1.7385377121166	0.024810520887464
-1.2340762153953	0.10901720602002
-0.73747372854539	0.28667550536283
-0.2453407083009	0.46224366960061
0.2453407083009	0.46224366960061
0.73747372854539	0.28667550536283
1.2340762153953	0.10901720602002
1.7385377121166	0.024810520887464
2.2549740020893	0.003243773342238
2.7888060584281	2.283386360164E-4
3.3478545673832	7.80255647853E-6
3.9447640401156	1.086069370769E-7
4.6036824495507	4.399340992273E-10
5.3874808900112	2.229393645534E-13];
w_SVDbasis = GH_int(:,2);

% Create data for this study
x = GH_int(:,1); x = pickpoints(-1,1,20);
N = length(x);
Neval = 100;
xeval = pickpoints(min(x),max(x),Neval);

% Choose an analytic function to study
yf = @(x) cos(x);yf = @(x) cos(3*pi*x);
y = yf(x);
yeval = yf(xeval);

% Choose a range of shape parameters
epvec = logspace(-1,1,20);
errvec = zeros(size(epvec));
errvecN = zeros(size(epvec));
errvecS = zeros(size(epvec));
errvecE = zeros(size(epvec));
errvecH = zeros(size(epvec));

k = 1;
for ep=epvec
    % Form the standard basis
    K = rbf(ep,DistanceMatrix(x,x));
    Kpos = K+1e-14*eye(N);
    Keval = rbf(ep,DistanceMatrix(xeval,x));
    warning('off','MATLAB:nearlySingularMatrix')
    errvec(k) = errcompute(Keval*(K\y),yeval);
    warning('on','MATLAB:nearlySingularMatrix')
    
    % Form the Newton basis
    chol_transform = chol(Kpos,'lower');
    Nmat = chol_transform;
    Neval = Keval/chol_transform';
    warning('off','MATLAB:nearlySingularMatrix')
    errvecN(k) = errcompute(Neval*(Nmat\y),yeval);
    warning('on','MATLAB:nearlySingularMatrix')
    
    % Form the SVD basis with weights 1
    % Note the application of the diagonal matrices is similar to HS-SVD
%     whalfvec = sqrt(w_SVDbasis);
%     whalfvec = ones(size(y));
%     [Q,Sigma2] = svd(bsxfun(@times,whalfvec,whalfvec').*K);
%     sigmavec = sqrt(diag(Sigma2));
%     Smat = Q.*bsxfun(@ldivide,whalfvec,sigmavec');
%     Seval = Keval*(Q.*bsxfun(@rdivide,whalfvec,sigmavec'));
%     warning('off','MATLAB:nearlySingularMatrix')
%     errvecS(k) = errcompute(Seval*(Smat\y),yeval);
%     warning('on','MATLAB:nearlySingularMatrix')
    
    % Form the eigenfunction basis with all N eigenfunctions
    gqr_alpha = 3;
    errvecE(k) = errcompute(gqr_eval(gqr_rsolve(x,y,ep,gqr_alpha,N),xeval),yeval);
    
    % Solve with the HS-SVD basis
    errvecH(k) = errcompute(gqr_eval(gqr_solve(x,y,ep,gqr_alpha),xeval),yeval);
    
    k = k + 1;
end

errpoly = errcompute(polyval(polyfit(x,y,N-1),xeval),yeval);

h = figure;
loglog(epvec,errvec,'--k','linewidth',2)
hold on
loglog(epvec,errvecN,'-sr','linewidth',2)
% loglog(epvec,errvecS,'-^m','linewidth',2)
loglog(epvec,errvecE,'-g','linewidth',2)
loglog(epvec,errvecH,'-+b','linewidth',2)
loglog(epvec,ones(size(epvec))*errpoly,':k','linewidth',2)
hold off
xlabel('$\varepsilon$','interpreter','latex')
ylabel('max norm error')
ylim([1e-8,1e4])
legend('standard','newton','eigenvalue','HS-SVD','polynomial','location','north')
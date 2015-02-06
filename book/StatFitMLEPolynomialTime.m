% StatFitMLEPolynomialTime
% Timing tests involving the StatFitMLEPolynomial example

% Load the data from the GaussQR web server
gqr_downloaddata('glacier_data.mat')
load glacier_data
[x,h] = rescale_data(datasites,heights,1);
N = size(x,1);

% Define some points at which to plot the solution
Nplot = 40;
xplot = pick2Dpoints(-1,1,Nplot);
X = reshape(xplot(:,1),Nplot,Nplot);
Y = reshape(xplot(:,2),Nplot,Nplot);

% Define an anisotropic missing Wendland kernel in dense form
rbf = @(r) (1-7*r.^2-81/4*r.^4).*sqrt(1-r.^2) - ...
           15/4*r.^4.*(6+r.^2).*log(r./(1+sqrt(1-r.^2)) + eps);

% Define a shape parameter vector
% We will fix this here and consider different polynomials
ep = [20 21];%ep = [5 4];

% Choose a variety of polynomial degrees to test
pd = [8;8];pd = [0;0];

% Define a function to evaluate a polynomial x(1).^a(1).*x(2).^a(2)
peval = @(x,a) prod(bsxfun(@power,x,a),2);
Pmat = @(x,arr) cell2mat(cellfun(@(a)peval(x,a),arr,'UniformOutput',0)');

% Number of runs to average over
nruns = 10;
time_factor = zeros(1,nruns);
time_mple = zeros(1,nruns);
time_ploteval = zeros(1,nruns);
for k=1:nruns
    % Prepare the Kriging matrix for use in choosing a polynomial
    tic
    K = DistanceMatrix(x,x,ep,rbf);
    [Lp,~,p] = chol(K,'lower','vector');
    logdetK = 2*sum(log(full(diag(Lp))));
    time_factor(k) = toc;
    
    tic
    alldegrees = num2cell((gqr_formMarr(pd+1)-1)',2);
    P = Pmat(x,alldegrees);
    LinvP_permuted = full(Lp\P(p,:));
    Linvh_permuted = Lp\h(p);
    PtKP = LinvP_permuted'*LinvP_permuted;
    beta = PtKP\(LinvP_permuted'*Linvh_permuted);
    c = zeros(N,1);
    c(p) = Lp'\(Linvh_permuted - LinvP_permuted*beta);
    mplevec = N*log((h-P*beta)'*c) + logdetK;
    time_mple(k) = toc;
    
    % Create a function to evaluate the universal kriging predictor
    tic
    sf = @(xe) DistanceMatrix(xe,x,ep,rbf)*c + Pmat(xe,alldegrees)*beta;
    splot = sf(xplot);
    time_ploteval(k) = toc;
end

fprintf('factoring = %g\nlikelihood = %g\nploteval = %g\n',mean(time_factor),mean(time_mple),mean(time_ploteval))
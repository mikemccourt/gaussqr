% StatFitMLEPolynomial
% This example explores the idea of using the MLE criterion to choose a
% polynomial term to accompany the kernel fit in universal kriging
% In this example, we study glacial data provided by Richard Franke
% Because of the size of this data (8345 data locations), we once again
% refer to the compactly supported Wendland kernels

% Load the data from the GaussQR web server
gqr_downloaddata('glacier_data.mat')
load glacier_data
[x,h] = rescale_data(datasites,heights,1);
N = size(x,1);
cind = convhull(x(:,1),x(:,2));

% Define some points at which to plot the solution
Nplot = 40;
xplot = pick2Dpoints(-1,1,Nplot);
X = reshape(xplot(:,1),Nplot,Nplot);
Y = reshape(xplot(:,2),Nplot,Nplot);

% Define an anisotropic missing Wendland kernel in dense form
rbf = @(r) (1-7*r.^2-81/4*r.^4).*sqrt(1-r.^2) - ...
           15/4*r.^4.*(6+r.^2).*log(r./(1+sqrt(1-r.^2)) + eps);
% Interesting side note: the following code produces complex results
% xplot = pick2Dpoints(-1,1,30);
% K = DistanceMatrix(xplot,z,[1 4],rbf);
% norm(imag(K),'fro')

% Define a shape parameter vector
% We will fix this here and consider different polynomials
ep = [5 4];ep = [8 11];ep = [16 14];ep = [20 21];

% Prepare the Kriging matrix for use in choosing a polynomial
K = DistanceMatrix(x,x,ep,rbf);return
[Lp,~,p] = chol(K,'lower','vector');
logdetK = 2*sum(log(full(diag(Lp))));

% Choose a variety of polynomial degrees to test
% pd2dvec = gqr_formMarr([6;6]) - 1;
pd1dvec = 0:8;
[P1,P2] = meshgrid(pd1dvec,pd1dvec);
pd2dvec = [P1(:)';P2(:)'];

% Define a function to evaluate a polynomial x(1).^a(1).*x(2).^a(2)
peval = @(x,a) prod(bsxfun(@power,x,a),2);
Pmat = @(x,arr) cell2mat(cellfun(@(a)peval(x,a),arr,'UniformOutput',0)');

% For each of the polynomials, compute the MPLE
k = 1;
mplevec = zeros(1,size(pd2dvec,2));
for pd=pd2dvec
    alldegrees = num2cell((gqr_formMarr(pd+1)-1)',2);
    P = Pmat(x,alldegrees);
    
	% Compute the polynomial mean term
    LinvP_permuted = full(Lp\P(p,:));
    Linvh_permuted = Lp\h(p);
    PtKP = LinvP_permuted'*LinvP_permuted;
    beta = PtKP\(LinvP_permuted'*Linvh_permuted);
    
    % Compute the kernel coefficients K\(h-P*beta)
    c = zeros(N,1);
    c(p) = Lp'\(Linvh_permuted - LinvP_permuted*beta);
    
    % Evaluate the likelihood
    mplevec(k) = N*log((h-P*beta)'*c) + logdetK;
    
    % Create a function to evaluate the universal kriging predictor
    sf = @(xe) DistanceMatrix(xe,x,ep,rbf)*c + Pmat(xe,alldegrees)*beta;
    
    k = k + 1;
end

% These plot commands only work for a 9x9 polynomial testing
h_mple = figure;
h_bar = bar3(pd1dvec,reshape(mplevec,size(P1)));
zlim([9e4,1.4e5])
xlim([0,10])
ylim([-1,9])
xlabel('2nd dimension degree')
ylabel('1st dimension degree')
zlabel('C_{MPLE}')
view([37.5 20])
h_submple = get(h_mple,'children');
set(h_submple,'xtick',[1 3 5 7 9])
set(h_submple,'ytick',[0 2 4 6 8])
set(h_submple,'xticklabel',{'0' '2' '4' '6' '8'})
arrayfun(@(n) set(h_bar(n),'facecolor',[.7 .7 .7]),1:9);

% Plot the final kriging prediction
h_pred = figure;
splot = sf(xplot);
S = reshape(splot,size(X));
S(not(inpolygon(X,Y,x(cind,1),x(cind,2)))) = NaN;
surf(X,Y,S,'edgealpha',.5)
zlim([1200 2100])
xlim([-1 1])
ylim([-1 1])
view([31 38])
colormap gray;C = colormap;colormap(flipud(C));
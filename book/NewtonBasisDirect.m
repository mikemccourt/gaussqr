% NewtonBasisDirect
% Computes and plots Newton basis functions for 2D RBF interpolation

% Consider the Gaussian RBF
rbf = @(e,r) exp(-(e*r).^2); ep = 3;

% Create some kernel centers and evaluation points
Ncenters = 5;
x = pick2Dpoints(0,1,Ncenters);
Neval = 40;
xeval = pick2Dpoints(0,1,Neval);

% Choose which Newton basis functions to plot
Nplot = [25 7 10];

% Evaluate the kernel matrices
DM_data = DistanceMatrix(x,x);
K = rbf(ep,DM_data);
DM_B = DistanceMatrix(xeval,x);
Keval = rbf(ep,DM_B);

% Find all Newton functions directly without iteration (works only for
% positive definite kernels)
N = chol(K,'lower');  % produces lower triangular matrix L
newtonbasisplot = Keval/N';

% Reshape the data for surface plotting
X = reshape(xeval(:,1),Neval,Neval);
Y = reshape(xeval(:,2),Neval,Neval);
Nmat = cellfun(@(nvec) reshape(nvec,Neval,Neval), ...
                      num2cell(newtonbasisplot,1),'UniformOutput',0);

% Create surface plots of the desired newton basis functions
h_surf = zeros(length(Nplot));
for k=1:length(Nplot)
    h_surf(k) = figure;
    surf(X,Y,Nmat{Nplot(k)},'FaceColor','interp','EdgeColor','none')
    colormap autumn; camlight; lighting gouraud
end
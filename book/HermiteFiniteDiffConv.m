% HermiteFiniteDiffConv
% This example studies the convergence of kernel-based finite difference
% methods for approximating derivatives with compact stencils.
% Specifically, we study the convergence for an increasing stencil size but
% fixed domain and show that the more points included in the stencil the
% better the approximate derivative is - up to a point.

% Choose to use the absolute max-norm error
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

% Choose a function to test
uf = @(x) exp(-2*x(:,1).*x(:,2));
ufL = @(x) 4*(x(:,1).^2+x(:,2).^2).*exp(-2*x(:,1).*x(:,2));

% Pick some points on which to provide data
% This should be relatively large
N = 100;
x = pick2Dpoints(0,1,N,'halt');
u = uf(x);

% Pick some points on which to evaluate the error
Neval = 25;
xeval = pick2Dpoints(.05,.95,Neval);
uLeval = ufL(xeval);

% Choose the C4 Matern kernel for computation
rbf = @(e,r) (1+(e*r)+(e*r).^2/3).*exp(-e*r);
rbfL = @(e,r) 1/3*e^2*exp(-e*r).*((e*r).^2-2*e*r-2);
ep = 1;

% Choose the various stencil sizes for convergence test
Nxvec = [5,12,25,50,80,130,200];

% Loop through the stencil sizes
% Compute and store the finite difference matrices
LFD = cell(size(Nxvec));
k = 1;
for Nx=Nxvec
    % Find the nearest neighbors to the evaluation points
    nearest = num2cell(knnsearch(x,xeval,'K',Nx),2);
    
    % Find the finite difference coefficients
    FDcell = cellfun(@(xe,xi)rbfL(ep,DistanceMatrix(xe,x(xi,:)))/...
                             rbf(ep,DistanceMatrix(x(xi,:),x(xi,:))),...
                      num2cell(xeval,2),nearest,'UniformOutput',0);

    % Form the vectors with the values for the sparse matrix
    FDvecs = cell2mat(cellfun(@(row,cols,vals)[row*ones(1,Nx);cols;vals],...
                      num2cell((1:Neval^2)',2),nearest,FDcell,'UniformOutput',0)');
    
	% Create the sparse matrix from those vectors
    LFD{k} = sparse(FDvecs(1,:),FDvecs(2,:),FDvecs(3,:),Neval^2,N^2);
    
    k = k + 1;
end

% Compute the error of the kernel FD Jacobian
errvec = cellfun(@(FDmat)errcompute(FDmat*u,uLeval),LFD);

h = figure;
loglog(Nxvec,errvec)
xlabel('stencil size')
ylabel('absolute max norm error')
% surf(reshape(xeval(:,1),Neval,Neval),reshape(xeval(:,2),Neval,Neval),reshape(uLeval,Neval,Neval))
% surf(reshape(xeval(:,1),Neval,Neval),reshape(xeval(:,2),Neval,Neval),reshape(A{1}*u,Neval,Neval))
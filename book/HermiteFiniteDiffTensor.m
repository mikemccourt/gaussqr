% HermiteFiniteDiffTensor
% This demonstrates that computing kernel based finite differences for
% nonradial kernels is possible

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

% Define tensor product kernels of the 1D C4 Materns
Kf = @(e,x,z) (1+e(1)*abs(bsxfun(@minus,x(:,1),z(:,1)'))+ ...
                (e(1)*abs(bsxfun(@minus,x(:,1),z(:,1)'))).^2/3).* ...
              exp(-e(1)*(abs(bsxfun(@minus,x(:,1),z(:,1)')))).* ...
              (1+e(2)*abs(bsxfun(@minus,x(:,2),z(:,2)'))+ ...
                (e(2)*abs(bsxfun(@minus,x(:,2),z(:,2)'))).^2/3).* ...
              exp(-e(2)*(abs(bsxfun(@minus,x(:,2),z(:,2)'))));
Kfxx = @(e,x,z) (e(1)^2/3).*(-1-e(1)*abs(bsxfun(@minus,x(:,1),z(:,1)'))+ ...
                               (e(1)*abs(bsxfun(@minus,x(:,1),z(:,1)'))).^2).* ...
                          exp(-e(1)*(abs(bsxfun(@minus,x(:,1),z(:,1)')))).* ...
                (1+e(2)*abs(bsxfun(@minus,x(:,2),z(:,2)'))+ ...
                  (e(2)*abs(bsxfun(@minus,x(:,2),z(:,2)'))).^2/3).* ...
                exp(-e(2)*(abs(bsxfun(@minus,x(:,2),z(:,2)'))));
Kfyy = @(e,x,z) (1+e(1)*abs(bsxfun(@minus,x(:,1),z(:,1)'))+ ...
                  (e(1)*abs(bsxfun(@minus,x(:,1),z(:,1)'))).^2/3).* ...
                 exp(-e(1)*(abs(bsxfun(@minus,x(:,1),z(:,1)')))).* ...
                (e(2)^2/3).*(-1-e(2)*abs(bsxfun(@minus,x(:,2),z(:,2)'))+ ...
                               (e(2)*abs(bsxfun(@minus,x(:,2),z(:,2)'))).^2).* ...
                          exp(-e(2)*(abs(bsxfun(@minus,x(:,2),z(:,2)'))));
KfL = @(e,x,z) Kfxx(e,x,z) + Kfyy(e,x,z);
epvec = [3,4];

% Choose the various stencil sizes for convergence test
Nxvec = [5,12,25,50,80,130,200];

% Loop through the stencil sizes
% Compute and store the finite difference matrices
A = cell(size(Nxvec));
k = 1;
for Nx=Nxvec
    % Find the nearest neighbors to the evaluation points
    nearest = num2cell(knnsearch(x,xeval,'K',Nx,'Distance','Chebychev'),2);
    
    % Find the finite difference coefficients
    FDcell = cellfun(@(xe,xi)KfL(epvec,xe,x(xi,:))/...
                             Kf(epvec,x(xi,:),x(xi,:)),...
                      num2cell(xeval,2),nearest,'UniformOutput',0);

    % Form the vectors with the values for the sparse matrix
    FDvecs = cell2mat(cellfun(@(row,cols,vals)[row*ones(1,Nx);cols;vals],...
                      num2cell((1:Neval^2)',2),nearest,FDcell,'UniformOutput',0)');
    
	% Create the sparse matrix from those vectors
    A{k} = sparse(FDvecs(1,:),FDvecs(2,:),FDvecs(3,:),Neval^2,N^2);
    
    k = k + 1;
end

% Compute the error of the kernel FD Jacobian
errvec = cellfun(@(FDmat)errcompute(FDmat*u,uLeval),A);

h = figure;
loglog(Nxvec,errvec)
xlabel('stencil size')
ylabel('absolute max norm error')
% NewtonBasisAdaptive
% Computes and plots Newton basis functions for 2D RBF interpolation

% Choose the Gaussians
rbf = @(e,r) exp(-(e*r).^2); ep = 5;

% Create kernel centers and evaluation points
Ncenters = 5;
x = pick2Dpoints(0,1,Ncenters);
Neval = 40;
xeval = pick2Dpoints(0,1,Neval);

% Create the kernel matrices
K = rbf(ep,DistanceMatrix(x,x));
Keval = rbf(ep,DistanceMatrix(xeval,x));

% This will be the matrix containing the Newton basis, columnwise,
% as functions on the evaluation grid.
newtonbasis = zeros(Ncenters^2,Ncenters^2);
newtonbasisplot = zeros(Neval^2,Ncenters^2);

% Create the first Newton basis function
z = diag(K);
[~,zind] = max(abs(z));
newtonbasis(:,1) = K(:,zind)/sqrt(z(zind));
newtonbasisplot(:,1) = Keval(:,zind)/sqrt(z(zind));

% Iterate through and compute one basis function at a time
w = zeros(Ncenters^2,1);
for iter=2:Ncenters
    % Determine the next kernel center to use
    % This is sort of a form of pivoting
    w = w + newtonbasis(:,iter-1).^2;
    [~,zind] = max(z-w);
    
    % Form the auxiliary terms for the Newton basis computation
    y = newtonbasis(zind,1:iter-1)';
    powerfuniter = sqrt(z(zind)-w(zind));
    
    % Evaluate the Newton basis on the kernel centers
    u = K(:,zind) - newtonbasis(:,1:iter-1)*y;
    newtonbasis(:,iter) = u/powerfuniter;
    
    % Evaluate the Newton basis on the plotting points
    u = Keval(:,zind) - newtonbasisplot(:,1:iter-1)*y;
    newtonbasisplot(:,iter) = u/powerfuniter;
end

% Plot the desired basis functions
figure
xe = reshape(xeval(:,1),Neval,Neval);
ye = reshape(xeval(:,2),Neval,Neval);
CFplot = surf(xe,ye,reshape(newtonbasisplot(:,ceil(Ncenters/2)),Neval,Neval));
set(CFplot,'FaceColor','interp','EdgeColor','none')
colormap autumn; view([145 45]); camlight; lighting gouraud
figure
CFplot = surf(xe,ye,reshape(newtonbasisplot(:,1),Neval,Neval));
set(CFplot,'FaceColor','interp','EdgeColor','none')
colormap autumn; view([145 45]); camlight; lighting gouraud
figure
CFplot = surf(xe,ye,reshape(newtonbasisplot(:,ceil(sqrt(Ncenters)/2)),Neval,Neval));
set(CFplot,'FaceColor','interp','EdgeColor','none')
colormap autumn; view([145 45]); camlight; lighting gouraud
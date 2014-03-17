% RBFFDConvergence.m
% This script plots some convergence behavior for RBF-FD

% Vector of points in the domain
Nvec = [50,50,50,50,50];

% Vector of # of points used to approximate the derivatives
kvec = [5,10,15,20,25];

% Vector of epsilon values used in the Gaussians
epvec = [1,1,1,1,1];

% Set to 1 if you want to plot point distribution
plot_points = 0;

% Loop through the requested values
errvec = zeros(size(Nvec));
for i=1:length(Nvec)
    errvec(i) = GQRFD_HelmholtzSolve(Nvec(i),kvec(i),epvec(i),plot_points);
end

% Plot the results
figure
semilogy(kvec,errvec)
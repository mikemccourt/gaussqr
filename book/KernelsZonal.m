% KernelsZonal
% This example demonstrates the computation of zonal kernels in Matlab
% through the computation of a dot-product term between data locations

% Create some points to compute on with the Matlab function sphere
% sphere returns points for plotting purposes, so the data must be reshaped
% in order to allow for evaluation, i.e., one point per row
Nplot = 105;
[X,Y,Z] = sphere(Nplot-1);
x = [X(:),Y(:),Z(:)];

% Choose a point at which to center the kernel
% Because we are using a color plot there is no easy way to plot more than
% one kernel, although computing with more than one kernel is possible
z = [1 0 2]/norm([1 0 2]);

% Define the inverse multiquadric zonal kernel
% This has a slightly different definition than the normal inverse
% multiquadric, and the shape parameter is bounded with g in [0,1].
% The dp term is not the distance r, but rather the pairwise dot product
% between the data locations x and z: this satisfies r^2 = 2*(1-x'*z).
% As before, each row represents one x location, each column is one z
zbf = @(g,dp) 1./sqrt(1+g^2-g*dp);
gam = .7;

% Evaluate the kernel
K = zbf(gam,x*z');

% Plot the kernel on a spherical surface
% This will require reshaping of the vector of kernel values into a matrix
% with the same dimensions as X, Y and Z
h = figure;
Kplot = reshape(K,Nplot,Nplot);
surf(X,Y,Z,Kplot,'edgecolor','none')

% Some extra stuff to help make the plot look nicer
axis square
view([1 1 1])
colorbar
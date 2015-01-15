% KernelsRadial2.m
% This is a demonstration of how to evaluate radial kernels in 2D
% The first example demonstrates isotropic kernels and the second
% demonstrates anisotropic kernels

% Choose some points at which to perform the evaluations
N = 50;
x = pick2Dpoints([-1,-1],[1,1],N);

% This surface plot requires reshaping of the data, which we do now rather
% than later
X = reshape(x(:,1),N,N);
Y = reshape(x(:,2),N,N);

% Choose some kernel centers to define the functions to be evaluated
z = [.6 .6;-.6 -.6;.6 -.6;-.6 .6;0 0];

% Compute the distance between the evaluation points and kernel centers
% DM = DistanceMatrix(x,z);
DM = pdist2(x,z);

% Define the radial kernel, here the inverse multiquadric
rbf_iso = @(e,r) 1./sqrt(1+(e*r).^2);

% Define shape parameter for the radial kernel
ep = 4;

% Kernel evaluation matrix, using pointwise evaluation of the distance
% matrix
K_iso = rbf_iso(ep,DM);

% Plotting the kernels on the domain
% This surface plot requires reshaping of the data
K1 = reshape(K_iso(:,1),N,N);
K2 = reshape(K_iso(:,2),N,N);
K3 = reshape(K_iso(:,3),N,N);
K4 = reshape(K_iso(:,4),N,N);
K5 = reshape(K_iso(:,5),N,N);
h_iso = figure;
hold on
surf(X,Y,K1,'edgecolor','none')
surf(X,Y,K2,'edgecolor','none')
surf(X,Y,K3,'edgecolor','none')
surf(X,Y,K4,'edgecolor','none')
surf(X,Y,K5,'edgecolor','none')
hold off
colormap jet
view([-.2 -1 1.3])
return

%%%%%%%%%%%
% Consider now the anisotropic kernel evaluated with the same centers

% Choose anisotropic shape parameters
epvec = [7,2];

% Compute the distance between the evaluation points and kernel centers
% Here we pass the shape parameter vector in
% The bug in Distance Matrix prevents me from running this as I should, so
% it's being run in this relatively dumb way right now.
DM_aniso = DistanceMatrix(x,z,epvec);

% Define the radial kernel, here the inverse multiquadric
rbf_aniso = @(r) 1./sqrt(1+r.^2);

% Kernel evaluation matrix, using pointwise evaluation of the distance
% matrix
K_aniso = rbf_aniso(DM_aniso);

% Plotting the kernels on the domain
K1 = reshape(K_aniso(:,1),N,N);
K2 = reshape(K_aniso(:,2),N,N);
K3 = reshape(K_aniso(:,3),N,N);
K4 = reshape(K_aniso(:,4),N,N);
K5 = reshape(K_aniso(:,5),N,N);
h_aniso = figure;
surf(X,Y,K1,'edgecolor','none')
hold on
surf(X,Y,K2,'edgecolor','none')
surf(X,Y,K3,'edgecolor','none')
surf(X,Y,K4,'edgecolor','none')
surf(X,Y,K5,'edgecolor','none')
hold off
colormap jet
view([-.2 -1 1.3])
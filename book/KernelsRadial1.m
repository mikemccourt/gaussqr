% KernelsRadial1.m
% This is a basic demonstration of how to evaluate radial kernels

% Choose some points at which to perform the evaluations
x = pickpoints(0,1,500);

% Choose some kernel centers to define the functions to be evaluated
z = [0;.4;.9];

% Compute the distance between the evaluation points and kernel centers
DM = DistanceMatrix(x,z);

% Define the radial kernel, here the C4 Matern
rbf = @(e,r) (1+e*r+(e*r).^2/3).*exp(-e*r);

% Define shape parameter for the radial kernel
ep = 5;

% Kernel evaluation matrix, using pointwise evaluation of the distance
% matrix
K = rbf(ep,DM);

% Plotting the kernels on the domain
plot(x,K,'linewidth',3)
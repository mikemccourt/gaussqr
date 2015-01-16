% KernelsTensor
% This example demonstrates how to compute tensor product kernels and how
% they are not the same as radial kernels in general
% We do this by comparing the radial C0 Matern kernel and the tensor
% product of 2 1D C0 Matern kernels
% This code relies heavily on the use of cell arrays to manage data from 1D
% while still being able to work in higher dimensions

% Define the C0 Matern kernel for 1D only
% Note that, although this kernel is radial, we implement it in the general
% x, z form to demonstrate how general kernels can be defined
% As always, each row is an x location, each column is a z location
rbf_1d = @(e,x,z) exp(-e*abs(bsxfun(@minus,x,z')));
epcell = {1.5,2.5};

% Choose some points in 1D to evaluate the tensor product kernel
% Store the data in cell arrays to allow varying N1d sizes
% cellfun performs the same operation on each term in an array
N1d = {101,99};
x1d = cellfun(@(x)pickpoints(-1,1,x),N1d,'UniformOutput',0);

% Choose z points at which to center this tensor product kernel
% This will also need to be stored as an array for later
Nz = 11;
z = pickpoints(0,1,Nz)*[1 6/5];
z1d = num2cell(z,1);

% Compute and store all the 1D kernel computations in a cell array
K1d = cellfun(@(e,x,z) rbf_1d(e,x,z),epcell,x1d,z1d,'UniformOutput',0);

% Form the full kernel matrix of size prod(cell2mat(N1d))-by-Nz
% In more than two dimensions this would likely require a loop invoking
% multiple calls to the kron function or structs or something
K = cell2mat(cellfun(@(Kx,Ky)kron(Kx,Ky), ...
    num2cell(K1d{1},1),num2cell(K1d{2},1),'UniformOutput',0));

% Define the C0 Matern radial kernel (same ep as before)
rbf = @(e,r) exp(-(e*r).^2);

% Pick a center to compare for plotting purposes
z = [-.6 -.5];

% Choose some 2D points to evaluate the radial kernel
N2d = 100;
x2d = pick2Dpoints([-1 -1],[1 1],[N2d N2d]);

% Evaluate the kernel
Kradial = rbf(ep,DistanceMatrix(x2d,z));

% Reshape the data so it can be surface plotted
[X1d,Y1d] = meshgrid(x1d{1},x1d{2});
X2d = reshape(x2d(:,1),N2d,N2d);
Y2d = reshape(x2d(:,2),N2d,N2d);
K_surf = reshape(K(:,6),N1d{2},N1d{1});
Kradial_surf = reshape(Kradial,N2d,N2d);

% Plot the tensor product and radial kernels on the same axes
surf(X2d,Y2d,Kradial_surf,'edgecolor','none')
hold on
surf(X1d,Y1d,K_surf,'edgecolor','none')
hold off
view([-1 1 1.5])
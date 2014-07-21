% ex22_gqr
% This is a basic example which demonstrates how RBF-FD works
% In essence, all that is required is a bunch of differentiation matrices
% which get compiled into one big FD operator matrix.
% In the simplest form, some neighboring points are used to approximate the
% derivative at a nearby point.
% Each point gets its own differentiation matrix, which is easy to do in a
% loop.  I'm working on how to do this in a vectorized structure.

% The problem we are trying to solve here is a Poisson equation
%       div * grad(u) = exp(x) + 2
%       u = exp(x) + x^2
interior_func = @(x) exp(x) + 2;
boundary_func = @(x) exp(x) + x.^2;
true_solution = @(x) exp(x) + x.^2;

% Choose the RBF that I want to use
% It'll be Gaussian so that I can compare to the GaussQR method
ep = .2;
rbf = @(e,r) exp(-(e*r).^2);
rbfxx = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);

% Choose the RBF-Direct or GaussQR solver
use_gaussqr = 1;

% Choose some points in the domain that are not evenly spaced, to emphasize
% that this does not require the standard FD structure
% N is just a guidance value that will be overwritten later
N = 40;
x_int = pickpoints(-1,1,N,'halton');
N_int = length(x_int);

% Add some points on the boundary, very easy in 1D
x_bc = [-1;1];
N_bc = length(x_bc);

% Redefine N as the total points in our problem
N = N_int + N_bc;

% Create the whole domain and associated access indices
x = [x_int;x_bc];
i_int = 1:N_int;
i_bc = N_int+1:N_int+N_bc;

% Loop through and figure out which points we want to use to approximate
% the differential operator
% There is likely a smarter way to do this, but this is the easiest
% After we find the nearest points, we compute the coefficients associated
% with approximating the derivative at that point
% Each step through this loop computes a 1 row differentiation matrix for
% evaluating the derivative at x(i) using the nearest points
n_closest = 6;
closest_indices = cell(N_int,1);
closest_coefs = cell(N_int,1);
for i=i_int
    this_x = x(i);
    [tmp,sorted_indices] = sort(abs(x-this_x));
    closest_indices{i} = sorted_indices(1:n_closest)';
    x_closest = x(sorted_indices(1:n_closest));
    
    % This step below performs the standard RBF solve
    % It just looks different because instead of solving Kc=y we are
    % solving it as y' = c'K  -->  c' = y'inv(K)
    % K is the interpolation matrix, sometimes called the Gram matrix
    % The right hand side is designed so that the differentiation matrix
    % annihilates our basis
    if use_gaussqr
        GQR = gqr_solveprep(0,x_closest,ep,1);
        psixx = gqr_phi(GQR,this_x,2)*[eye(n_closest);GQR.Rbar];
        Psi = GQR.stored_phi1 + GQR.stored_phi2*GQR.Rbar;
        Diff_Mat_this_x = psixx/Psi;
    else
        DistMat = DistanceMatrix(x_closest,x_closest);
        K = rbf(ep,DistMat);
        DiffOpMat = DistanceMatrix(this_x,x_closest);
        kxx = rbfxx(ep,DiffOpMat);
        Diff_Mat_this_x = kxx/K;
    end
    closest_coefs{i} = Diff_Mat_this_x;
end

% Define the matrix that we need to invert for this problem: Ax=b
A = zeros(N);

% Now that we have the various differentiation matrices, we can compose our
% full differential operator
% We obviously could have put all this together above, but I want to
% isolate these two steps to demonstrate explicitly what's happening
% The loop below fills each row with the differentiation matrix
% (coefficients) associated with that point
for i=i_int
    A(i,closest_indices{i}) = closest_coefs{i};
end

% We also must define the boundary condition operator
% In this case, with Dirichlet boundary conditions, it's just the identity
% Rows from the boundary condition operator are placed in A to enforce the
% boundary conditions
BC_operator = eye(N);
for i=i_bc
    A(i,:) = BC_operator(i,:);
end

% The RHS for this problem is defined by the functions earlier
b = zeros(N,1);
b(i_int) = interior_func(x(i_int));
b(i_bc) = boundary_func(x(i_bc));

% Solve the system and plot the error
u = A\b;
[x_plot,plot_ind] = sort(x);
u_plot = u(plot_ind);
err_plot = abs(u_plot-true_solution(x_plot));
subplot(1,3,1)
semilogy(x_plot,err_plot,'linewidth',3)
ylabel('Absolute Error')
title(sprintf('N=%d, k=%d, ep=%g',N,n_closest,ep))

% Plot the sparsity pattern for this matrix
% Note that the sparsity pattery here is not the optimal one, it's just the
% one that's generated based on the very straightforward code we have
% written here.  A better ordering of the points would yield a better
% bandwidth, as demonstrated below.
subplot(1,3,2)
spy(A)
title('Simple Ordering')
subplot(1,3,3)
r = symrcm(A);
spy(A(r,r))
title('RCM Ordering')
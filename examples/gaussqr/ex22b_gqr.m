% ex22b_gqr
% This is a 2D example which demonstrates how RBF-FD works
% In essence, all that is required is a bunch of differentiation matrices
% which get compiled into one big FD operator matrix.
% In the simplest form, some neighboring points are used to approximate the
% derivative at a nearby point.
% Each point gets its own differentiation matrix, which is easy to do in a
% loop.  I'm working on how to do this in a vectorized structure.

% The problem we are trying to solve here is a Helmholtz equation
%     Lap(u) - lambda^2*u = f   -on- interior
%     u = g                     -on- boundary
%   solution: u(x,y) = sin(x^2+y)
%   domain: L shaped region (-1<x<0 & -1<y<1)+(-1<x<1 & -1<y<0)
% On the y=1 boundary we will use Neumann BC, Dirichlet everywhere else
% lambda is the wavenumber (I think) for this problem
lambda = 3;
interior_func = @(x,y) 2*cos(x.^2+y)-4*x.^2.*sin(x.^2+y)-sin(x.^2+y)-lambda^2*sin(x.^2+y);
dirichlet_bc  = @(x,y) sin(x.^2+y); % Also the true solution
neumann_bc    = @(x,y) cos(x.^2+y);

% Choose the RBF that I want to use
% It'll be Gaussian so that I can compare to the GaussQR method
ep = 1;
rbf = @(e,r) exp(-(e*r).^2);
rbfy = @(e,r,dy) -2*e^2*dy.*exp(-(e*r).^2);
rbfL = @(e,r) 4*e^2*((e*r).^2-1).*exp(-(e*r).^2);

% Choose the RBF-Direct or GaussQR solver
% NOTE: Runs slowly, but it seems to work fine
use_gaussqr = 0;

% Choose some points in the domain that are not evenly spaced, to emphasize
% that this does not require the standard FD structure
% N is just a 1D surrogate for however many points we'll actually be
% working with, which will be closer to 3/4N^2-4N
N = 30;
% Find the whole square and dump the top-right quadrant
x_halt = pick2Dpoints([-1,-1],[1,1],N,'halton');
x_int = x_halt(x_halt(:,1)<0 | x_halt(:,2)<0,:);
N_int = length(x_int);

% Add some points on the boundary, a little harder in 2D
% Distinguish between the Dirichlet and Neumann boundaries
x_m11 = pickpoints(-1,1,N);
x_01 = pickpoints(0,1,N/2);
x_bcN = [x_m11(2:end-1),-ones(N-2,1)];
x_bcD = [[-ones(N,1),x_m11]; ...
         [-x_01(2:end-1),ones(N/2-2,1)]; ...
         [zeros(N/2,1),x_01]; ...
         [x_01(2:end-1),zeros(N/2-2,1)]; ...
         [ones(N/2,1),-x_01]];
N_bcN = length(x_bcN);
N_bcD = length(x_bcD);

% Redefine N as the total points in our problem
N = N_int + N_bcN + N_bcD;

% Create the whole domain and associated access indices
x = [x_int;x_bcN;x_bcD];
i_int = 1:N_int;
i_bcN = N_int+1:N_int+N_bcN;
i_bcD = N_int+N_bcN+1:N_int+N_bcN+N_bcD;

% Loop through and figure out which points we want to use to approximate
% the differential operator
% There is likely a smarter way to do this, but this is the easiest
% After we find the nearest points, we compute the coefficients associated
% with approximating the derivative at that point
% Each step through this loop computes a 1 row differentiation matrix for
% evaluating the derivative at x(i) using the nearest points
% Note that this loop includes both the finite difference computation for
% the interior Helmholtz operator and the Neumann boundary conditions
n_closest = 15;
closest_indices = cell(N_int+N_bcD,1);
closest_coefs = cell(N_int+N_bcD,1);
for i=[i_int,i_bcN]
    this_x = x(i,:);
    [tmp,sorted_indices] = sort(sum((x-repmat(this_x,N,1)).^2,2));
    closest_indices{i} = sorted_indices(1:n_closest)';
    x_closest = x(closest_indices{i},:);
    
    % This step below performs the standard RBF solve
    % It just looks different because instead of solving Kc=y we are
    % solving it as y' = c'K  -->  c' = y'inv(K)
    % K is the interpolation matrix, sometimes called the Gram matrix
    % The right hand side is designed so that the differentiation matrix
    % annihilates our basis
    % Also note that the differential operator is Laplacian - lambda^2*I
    if use_gaussqr
        GQR = gqr_solveprep(0,x_closest,ep,1);
        InterpMat = GQR.stored_phi1 + GQR.stored_phi2*GQR.Rbar;
        if ismember(i,i_int)
            DiffOp = (gqr_phi(GQR,this_x,[2,0])+gqr_phi(GQR,this_x,[0,2]))*[eye(n_closest);GQR.Rbar];
        else
            DiffOp = gqr_phi(GQR,this_x,[0,1])*[eye(n_closest);GQR.Rbar];
        end
    else
        DistMat = DistanceMatrix(x_closest,x_closest);
        InterpMat = rbf(ep,DistMat);
        DiffOpDistMat = DistanceMatrix(this_x,x_closest);
        if ismember(i,i_int)
            DiffOp = rbfL(ep,DiffOpDistMat) - lambda^2*rbf(ep,DiffOpDistMat);
        else
            yDiffMat = DifferenceMatrix(this_x(:,2),x_closest(:,2));
            DiffOp = rbfy(ep,DiffOpDistMat,yDiffMat);
        end
    end
    closest_coefs{i} = DiffOp/InterpMat;
end

% Define the matrix that we need to invert for this problem: Ax=b
A = zeros(N);

% Now that we have the various differentiation matrices, we can compose our
% full differential operator
% We obviously could have put all this together above, but I want to
% isolate these two steps to demonstrate explicitly what's happening
% The loop below fills each row with the differentiation matrix
% (coefficients) associated with that point
% Note that the coefficients we computed only evaluate the Laplacian
% We also need to apply the -lambda^2u
for i=i_int
    A(i,closest_indices{i}) = closest_coefs{i};
end

% We also must define the boundary condition operator
% In this case, with Dirichlet boundary conditions, it's just the identity
% Rows from the boundary condition operator are placed in A to enforce the
% boundary conditions
% We have already found the necessary Neumann operator above, and its
% contributions are included in the first loop below
for i=i_bcN
    A(i,closest_indices{i}) = closest_coefs{i};
end
bcD_operator = eye(N);
for i=i_bcD
    A(i,:) = bcD_operator(i,:);
end

% The RHS for this problem is defined by the functions earlier
% Note that to make the function definitions cleaner we explicitly separate
% the x and y components here and pass them separately
b = zeros(N,1);
b(i_int) = interior_func(x(i_int,1),x(i_int,2));
b(i_bcN) = neumann_bc(x(i_bcN,1),x(i_bcN,2));
b(i_bcD) = dirichlet_bc(x(i_bcD,1),x(i_bcD,2));

% Solve the system
u = A\b;

% Plot the pointwise error
% Note that we are going to plot the error on a log scale by first taking
% the log and interpolating that
subplot(1,3,1)
err_solve = log10(abs(u - dirichlet_bc(x(:,1),x(:,2))));
err_plot = TriScatteredInterp(x(:,1),x(:,2),err_solve);
x_plot = pickpoints(-1,1,60);
[X,Y] = meshgrid(x_plot,x_plot);
E = err_plot(X,Y);
surf(X,Y,E)
% u_plot = TriScatteredInterp(x(:,1),x(:,2),dirichlet_bc(x(:,1),x(:,2)));
% u_plot = TriScatteredInterp(x(:,1),x(:,2),u);
% U = u_plot(X,Y);
% surf(X,Y,U)
% zlim([-1,1])
plot(x_bcN(:,1),x_bcN(:,2),'or',x_bcD(:,1),x_bcD(:,2),'bx',x_int(:,1),x_int(:,2),'m.')
zlabel('log_{10} Absolute Error')
E_err = errcompute(u,dirichlet_bc(x(:,1),x(:,2)));
title(sprintf('N=%d, k=%d, ep=%g, err=%g',N,n_closest,ep,E_err))

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
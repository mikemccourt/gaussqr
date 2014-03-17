function [err,N,timevec] = GQRFD_HelmholtzSolve(N,k,ep,plot_points)
% This function computes the error to the solution of a Helmholtz problem
% we have specified:
%     Lap(u) - lambda^2*u = f   -on- interior
%     u = g                     -on- boundary
%   solution: u(x,y) = sin(x^2+y)
%   domain: L shaped region (-1<x<0 & -1<y<1)+(-1<x<1 & -1<y<0)
% On the y=1 boundary we will use Neumann BC, Dirichlet everywhere else
% lambda is the wavenumber (I think) for this problem
% The solution is computed with Gaussian generated finite differences
%
% function err = GQRFD_HelmholtzSolve(N,k,ep,plot_points)
%   Inputs : N - 1D surrogate for number of points (3/4N^2+4N)
%            k - Number of points to use for finite differences
%            ep - Gaussian shape parameter
%            plot_points - (optional, default=0) Plot the point distribution used
%  Outputs : err - The RMS error of the solution
%            N - The actual number of points that appear in the domain
%            timevec - The times that this program takes
%                      timevec(1) - time to find A values
%                      timevec(2) - time to construct and solve A\b

if nargin<4
    plot_points = 0;
end
timevec = zeros(1,2);

% This is the PDE description
lambda = 3;
interior_func = @(x,y) 2*cos(x.^2+y)-4*x.^2.*sin(x.^2+y)-sin(x.^2+y)-lambda^2*sin(x.^2+y);
dirichlet_bc  = @(x,y) sin(x.^2+y); % Also the true solution
neumann_bc    = @(x,y) cos(x.^2+y);

% Set up the Gaussian RBF and shape parameter
rbf = @(e,r) exp(-(e*r).^2);
rbfy = @(e,r,dy) -2*e^2*dy.*exp(-(e*r).^2);
rbfL = @(e,r) 4*e^2*((e*r).^2-1).*exp(-(e*r).^2);

% Choose some points in the domain that are Halton spaced
% Find the whole square and dump the top-right quadrant
point_generator = haltonset(2,'Skip',1);
x_halt = 2*net(point_generator,N^2) - 1;
x_int = x_halt(x_halt(:,1)<0 | x_halt(:,2)<0,:);
N_int = size(x_int,1);

% Add some points on the boundary, a little harder in 2D
% Distinguish between the Dirichlet and Neumann boundaries
x_m11 = linspace(-1,1,N)';
x_01 = linspace(0,1,N/2)';
x_bcN = [x_m11(2:end-1),-ones(N-2,1)];
x_bcD = [[-ones(N,1),x_m11]; ...
         [-x_01(2:end-1),ones(N/2-2,1)]; ...
         [zeros(N/2,1),x_01]; ...
         [x_01(2:end-1),zeros(N/2-2,1)]; ...
         [ones(N/2,1),-x_01]];
N_bcN = size(x_bcN,1);
N_bcD = size(x_bcD,1);

% Redefine N as the total points in our problem
N = N_int + N_bcN + N_bcD;

% Create the whole domain and associated access indices
x = [x_int;x_bcN;x_bcD];
i_int = 1:N_int;
i_bcN = N_int+1:N_int+N_bcN;
i_bcD = N_int+N_bcN+1:N_int+N_bcN+N_bcD;

% Prepare storage for the matrix
tic
rowvec = zeros((N_int + N_bcN)*k + N_bcD,1);
colvec = zeros((N_int + N_bcN)*k + N_bcD,1);
valvec = zeros((N_int + N_bcN)*k + N_bcD,1);
timevec(2) = toc;

% Loop through and figure out which points we want to use to approximate
% the differential operator:
% Either Laplacian - lambda^2*Identity or d/dy on the boundary
% We will store the matrix as a sparse matrix
tic
j = 1; % Index of vector location for sparse matrix
for i=[i_int,i_bcN]
    this_x = x(i,:);
    [~,sorted_indices] = sort(sum((x-repmat(this_x,N,1)).^2,2));
    closest_indices = sorted_indices(1:k)';
    x_closest = x(closest_indices,:);
    
    DistMat = DistanceMatrix(x_closest,x_closest);
    InterpMat = rbf(ep,DistMat);
    DiffOpDistMat = DistanceMatrix(this_x,x_closest);
    if ismember(i,i_int)
        DiffOp = rbfL(ep,DiffOpDistMat) - lambda^2*rbf(ep,DiffOpDistMat);
    else
        yDiffMat = DifferenceMatrix(this_x(:,2),x_closest(:,2));
        DiffOp = rbfy(ep,DiffOpDistMat,yDiffMat);
    end
    closest_coefs = DiffOp/InterpMat;
    
    for m=1:k
        rowvec(j) = i;
        colvec(j) = closest_indices(m);
        valvec(j) = closest_coefs(m);
        j = j + 1;
    end
end

% We still need to add the Dirichlet BC rows: just the identity
for i=i_bcD
    rowvec(j) = i;
    colvec(j) = i;
    valvec(j) = 1;
    j = j + 1;
end
timevec(1) = toc;

% Define the matrix that we need to invert for this problem: Ax=b
tic
A = sparse(rowvec,colvec,valvec);

% The RHS for this problem is defined by the functions earlier
b = zeros(N,1);
b(i_int) = interior_func(x(i_int,1),x(i_int,2));
b(i_bcN) = neumann_bc(x(i_bcN,1),x(i_bcN,2));
b(i_bcD) = dirichlet_bc(x(i_bcD,1),x(i_bcD,2));

% Solve the system and compute the RMS error
u = A\b;
err = norm(u-dirichlet_bc(x(:,1),x(:,2)))/sqrt(N);
timevec(2) = timevec(2) + toc;

% Plot the error if requested
if plot_points
    figure
    subplot(1,3,1)
    plot(x_bcN(:,1),x_bcN(:,2),'or',x_bcD(:,1),x_bcD(:,2),'bx',x_int(:,1),x_int(:,2),'m.')
    title(sprintf('N=%d, k=%d, ep=%g, err=%g',N,k,ep,err))
    subplot(1,3,2)
    spy(A)
    title('Simple Ordering')
    subplot(1,3,3)
    r = symrcm(A);
    spy(A(r,r))
    title('RCM Ordering')
end
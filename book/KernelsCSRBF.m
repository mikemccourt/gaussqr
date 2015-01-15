% KernelsCSRBF.m
% This is a brief example of computing with sparse distance matrics and the
% compactly supported kernels that exist as a result.

% These functions are used to manipulate the sparse matrices without
% causing them to be dense matrices
% sppc - sparse plus constant (although the arguments are in reverse order)
% spd  - sparse divided by sparse
sppc = @(const,Mat) spfun(@(x)const + x,Mat);
spd  = @(Mat1,Mat2) Mat1.*spfun(@(x)1./x,Mat2);

% Definition of the compactly supported kernel that can be applied to a
% sparse matrix without changing the structure
%
% Note that the shape parameter must be handled within the DistanceMatrix
% computation, thus the radial kernel is just a single argument
%
% Also note that, because the distance matrix will be sparse, the spfun
% computation in needed to avoid log(0) terms from occurring and to avoid
% 1+sqrt(1+r.^2) from introducing nonzero values throughout the matrix
%
% The term in the log is actually
%              r/(1+sqrt(1+eps-r^2))
% To avoid zeros being squeezed out of the matrix when r = 1 because when
% Matlab stores values in a matrix it dumps zero-valued nonzeros
%
% Another eps term appears added to the log argument to prevent a log(0)
% when r=0, which will get balanced out by the r^2 term
rbf = @(r) sppc(1,2*r.^2).*sqrt(sppc(1,-r.^2)) + 3*r.^2.* ...
           spfun(@log,sppc(eps,spd(r,sppc(1+eps,sqrt(sppc(1,-r.^2))))));
       
% Definition of the compactly supported kernel that can be applied within
% the DistanceMatrix function as though it were dense
% to avoid the weird sparse matrix operations used in rbf
%
% This is easier to code and implement, but is more costly if a single
% distance matrix can be formed and then evaluated with several different
% kernels for multiple tests
%
% The eps added to the log argument serves the same purpose as before;
% avoiding a log(0) term when r=0
rbf_in_DM = @(r) (1+2*r.^2).*sqrt(1-r.^2) + ...
                 3*r.^2.*log(r./(1+sqrt(1-r.^2))+eps);
             
% Choose a shape parameter or shape parameters
ep = 4;
epvec = [5 1];
             
% Choose some points in 2D at which to sample the kernels
N = 70;
x = pick2Dpoints([-1,-1],[1,1],N);

% This surface plot requires reshaping of the data, which we do now rather
% than later
X = reshape(x(:,1),N,N);
Y = reshape(x(:,2),N,N);

% Choose some kernel centers to define the functions to be evaluated
% z_in_DM will be evaluated by passing rbf_in_DM
z_in_DM = [-.6 .6;-.6 -.6;0 0];
z = [.6 -.6;.6 .6];

% Evaluate the sparse distance matrix, the kernels and plot them
% For this evaluation we will use an isotropic kernel
DM_iso = DistanceMatrix(x,z,ep,1);
K_iso = rbf(DM_iso);

% Plot the kernels we have computed thus far
K1 = reshape(K_iso(:,1),N,N);
K2 = reshape(K_iso(:,2),N,N);
hold on
surf(X,Y,K1,'edgecolor','none')
surf(X,Y,K2,'edgecolor','none')

% Compute some anisotropic kernels by passing the rbf function handle to
% the distance matrix function
K_aniso = DistanceMatrix(x,z_in_DM,epvec,rbf_in_DM);

% Plot these new results on the same axes
K3 = reshape(K_aniso(:,1),N,N);
K4 = reshape(K_aniso(:,2),N,N);
K5 = reshape(K_aniso(:,3),N,N);
surf(X,Y,K3,'edgecolor','none')
surf(X,Y,K4,'edgecolor','none')
surf(X,Y,K5,'edgecolor','none')

hold off
colormap jet
view([-.2 -1 1.3])
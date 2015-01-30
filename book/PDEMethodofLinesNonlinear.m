% PDEMethodofLinesNonlinear
% This solves a 1D heat equation with collocation in the spatial dimension
% and the method of lines in the time dimension.  The PDE is 
%       u_t = .1*u_xx + u - u^3,  u(-1) = -1,  u(1) = 1
%                                  u(x,0) = x + sin(pi*x)
% which has been adapted from [Trefethen (2000)]

% Create the initial, boundary and interior functions
fic = @(x) x - 3*sin(pi*x);
fbc = @(x) x;
fint = @(x) zeros(size(x,1),1);

% Choose to use the the C4 Matern kernels here
rbf = @(e,r) exp(-e*r).*(1+e*r+(e*r).^2/3);
rbfxx = @(e,r) e^2*exp(-e*r).*((e*r).^2-e*r-1)/3;
ep = 3;

% Choose the time range for this problem
tspan = [0,3];

% Choose some number of discretization points in the interior
% Create the boundary points and the full collocation vector
Nx = 20;
xint = pickpoints(-1,1,Nx,'halt');
xbc = [-1;1];
x = [xint;xbc];

% Create the necessary distance and collocation matrices
DMint = DistanceMatrix(xint,x);
DMbc = DistanceMatrix(xbc,x);
Kint = rbf(ep,DMint);
Kxxint = rbfxx(ep,DMint);
Kbc = rbf(ep,DMbc);
K = [Kint;Kbc];

% Form inhomogeneous term associated with the boundary
% Note that it is negative to enforce the BC and not -BC below
fgvec = [fint(xint);-fbc(xbc)];

% Define the ODE system function u' = f(t,u)
odefun = @(t,c) [.1*Kxxint+Kint;Kbc]*c - ([Kint;zeros(size(Kbc))]*c).^3 + fgvec;

% Define the Mass matrix, explaining which rows are boundary rows
Mass = [Kint;zeros(size(Kbc))];
odeopts = odeset('Mass',Mass,'MassSingular','yes');

% Create the initial conditions by interpolating the initial data
cinit = K\fic(x);

% Solve the ODE
[tvec,C] = ode15s(odefun,tspan,cinit,odeopts);

% Plot the solution
xeval = pickpoints(-1,1,30);
Keval = rbf(ep,DistanceMatrix(xeval,x));
U = Keval*C';
surf(xeval,tvec',U','edgecolor','none')

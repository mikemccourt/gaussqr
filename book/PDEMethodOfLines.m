% PDEMethodofLines
% This example solves a heat equation using the HS-SVD basis collocation in
% space and the method of lines in time.
%           u_t = .06*u_xx + exp(sin(3t)),     
%           u(0,t) = 0,    u(1,t) = 1,        u(x,0) = (1-x)^2

% Choose to use a kernel, in this case a C4 Wendland kernel
rbf = @(e,r) (35*(e*r).^2+18*e*r+3).*max(1-e*r,0).^6;
rbfxx = @(e,r) 56*e^2*(35*(e*r).^2-4*e*r-1).*max(1-e*r,0).^4;
ep = 3;

% Create the initial, boundary and interior functions
fic = @(x) (1-x).^2;
fbc = @(x) 1-x;
fint = @(x,t) ones(size(x,1),1)*exp(sin(3*t));

% Choose the time range for this problem
tspan = [0,10];

% Choose some number of discretization points in the interior
% Create the boundary points and the full collocation vector
Nx = 20;
xint = pickpoints(0,1,Nx,'halt');
xbc = [0;1];
x = [xint;xbc];

% Create the necessary distance and collocation matrices
DMint = DistanceMatrix(xint,x);
DMbc = DistanceMatrix(xbc,x);
Kint = rbf(ep,DMint);
Kxxint = rbfxx(ep,DMint);
Kbc = rbf(ep,DMbc);
K = [Kint;Kbc];

% Define the ODE system function u' = f(t,u)
odefun = @(t,c) [.06*Kxxint;Kbc]*c + [fint(xint,t);-fbc(xbc)];

% Define the Mass matrix, explaining which rows are boundary rows
Mass = [Kint;zeros(size(Kbc))];
odeopts = odeset('Mass',Mass,'MassSingular','yes','MStateDependence','none');

% Create the initial conditions by interpolating the initial data
cinit = K\fic(x);

% Solve the ODE
[teval,C] = ode15s(odefun,tspan,cinit,odeopts);

% Evaluate the solution
xeval = pickpoints(0,1,30);
Keval = rbf(ep,DistanceMatrix(xeval,x));
U = Keval*C';

% Plot the solution
h = figure;
surf(xeval,teval',U','edgecolor','none')
xlabel('x'),ylabel('t'),zlabel('u'),zlim([0,4]),view([-.3 -1 1])
% PDEMFSDemo
% MFS solver for a Helmholtz problem on a L-shaped domain with mixed
% boundary conditions
% Lap(u)+nu^2*u = 0                on interior, nu^2 = 2*pi^2
% d/dn(u) = -pi*sin(pi*(x+y))    on y = 2
% u = cos(pi*(x+y))              on other boundary

nu = sqrt(2)*pi;
% Define the kernel, here the fundamental solution
Gf = @(r) 1i/4*besselh(0,2,nu*r);
Gyf = @(r,dy) -1i*nu/4*besselh(1,2,nu*r).*dy./r;

% Set up the functions for the boundary conditions
fdir = @(x) cos(pi*(x(:,1)+x(:,2)));
fneu = @(x) -pi*sin(pi*(x(:,1)+x(:,2)));

% Define the number of points to use on each boundary
N = 15;

% Construct square in 2D from 2D tensor grid
xbox = pick2Dpoints(0,2,N);
xsquare = xbox(any(xbox==0 | xbox==2,2),:);
% Create L shape from square shape
xneedflipped = all(xsquare<=1,2);
xL = xsquare;  xL(xneedflipped,:) = bsxfun(@minus,[1 1],xsquare(xneedflipped,:));
% Find the Neumann and Dirichlet boundaries
xneu = xsquare(xL(:,2)==2,:);
xdir = setdiff(xL,xneu,'rows');

% Choose Green's kernel source points (the centers)
theta = pickpoints(0,2*pi,4*N);
z = bsxfun(@plus,[1 1],2*[cos(theta),sin(theta)]);

% Define the distance matrices
DMdir = DistanceMatrix(xdir,z);
DMneu = DistanceMatrix(xneu,z);
DiffMyneu = DifferenceMatrix(xneu(:,2),z(:,2));

% Define the kernel matrices
Gdir = Gf(DMdir);
Gneu = Gyf(DMneu,DiffMyneu);
A = [Gdir;Gneu];

% Define the right hand side
rhsdir = fdir(xdir);
rhsneu = fneu(xneu);
rhs = [rhsdir;rhsneu];

% Solve the system
coef = A\rhs;

% Evaluate the solution, and NaN points outside the domain
Neval = 21;
xeval = pick2Dpoints(0,2,Neval);
Geval = Gf(DistanceMatrix(xeval,z));
ueval = real(Geval*coef);
ueval(all(xeval<1,2)) = NaN;

% Plot the results on a surface plot
X = reshape(xeval(:,1),Neval,Neval);
Y = reshape(xeval(:,2),Neval,Neval);
U = reshape(ueval,Neval,Neval);
surf(X,Y,U),view([-.3 -2.1 7])

h_pts = figure;
plot(xdir(:,1),xdir(:,2),'or','linewidth',2)
hold on
plot(xneu(:,1),xneu(:,2),'+m','linewidth',2)
plot(z(:,1),z(:,2),'xb','linewidth',2)
hold off
legend('Dirichlet BC','Neumann BC','kernel centers','location','southwest')
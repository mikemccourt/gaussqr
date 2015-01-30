% PDEMPSDemo
% MPS solver for a Helmholtz problem on a L-shaped domain with mixed
% boundary conditions
% Lap(u)-nu^2*u = -exp(2*x+2*y)     on interior, nu = 3
% u = exp(2*x+2*y)                   on boundary
% The domain is defined by r(t) = .8 + .1*(sin(6*t)+sin(3*t))

nu = 3;
% Define the MFS kernel, here the fundamental solution
Gf = @(r) 1/(2*pi)*besselk(0,nu*r);

% Define the MPS kernel, in this case a C4 Wendland kernel
rbf = @(e,r) (35*(e*r).^2+18*e*r+3).*max(1-e*r,0).^6;
rbfL = @(e,r) 56*e^2*(40*(e*r).^2-8*e*r-2).*max(1-e*r,0).^4;
ep = .1;

% Set up the functions for the boundary conditions
fbc = @(x) exp(2*x(:,1)+2*x(:,2));
fint = @(x) -exp(2*x(:,1)+2*x(:,2));

% Define the number of points to use in each region
Nbc = 50;  Nint = 200;

% Choose collocation points on the boundary
theta = pickpoints(0,2*pi,Nbc);
bcdef = @(t) .8 + .1*(sin(6*t)+sin(3*t));
xbc = [bcdef(theta).*cos(theta),bcdef(theta).*sin(theta)];

% Choose MFS centers
zbc = 2*[cos(theta),sin(theta)];

% Choose collocation points/centers in interior
% Create a lot of Halton points and keep the first Nint
xtest = pick2Dpoints(-1,1,sqrt(3*Nint),'halt');
xallin = xtest(inpolygon(xtest(:,1),xtest(:,2),xbc(:,1),xbc(:,2)),:);
xint = xallin(1:Nint,:);

% Define the distance matrices
DMint = DistanceMatrix(xint,xint);
DMintformfs = DistanceMatrix(xbc,xint);
DMbc = DistanceMatrix(xbc,zbc);

% Define the MPS matrix and solve the MPS problem
K = rbf(ep,DMint);
KL = rbfL(ep,DMint);
A = KL - nu^2*K;
mpscoef = A\fint(xint);

% Define the MFS system and solve the MFS problem
G = Gf(DMbc);
Kformfs = rbf(ep,DMintformfs);
mfsrhs = fbc(xbc) - Kformfs*mpscoef;
mfscoef = G\mfsrhs;

% Define a function to evaluate the combined solution
usol = @(x)      Gf(DistanceMatrix(x,zbc))*mfscoef + ...
            rbf(ep,DistanceMatrix(x,xint))*mpscoef;

% Evaluate the solution, and NaN points outside the domain
Neval = 55;
xeval = pick2Dpoints(-1,1,Neval);
ueval = usol(xeval);
uerr = abs(usol(xeval) - fbc(xeval));
ueval(not(inpolygon(xeval(:,1),xeval(:,2),xbc(:,1),xbc(:,2)))) = NaN;

% Plot the results on a surface plot
h_surf = figure;
X = reshape(xeval(:,1),Neval,Neval);
Y = reshape(xeval(:,2),Neval,Neval);
U = reshape(ueval,Neval,Neval);
Uerr = reshape(uerr.*not(isnan(ueval))+ 1e-5*isnan(ueval),Neval,Neval);
surf(X,Y,U,log10(Uerr)),view([-5 47])
h_color = colorbar('fontsize',14);
set(h_color,'ticks',[-8,-6,-4])
set(h_color,'ticklabels',{'10^{-8}','10^{-6}','10^{-4}'})

h_pts = figure;
plot(xbc(:,1),xbc(:,2),'or')
hold on
plot(zbc(:,1),zbc(:,2),'xb')
plot(xin(:,1),xin(:,2),'.k','markersize',11)
hold off
legend('MFS collocation','MFS centers','MPS collocation','location','northeast')
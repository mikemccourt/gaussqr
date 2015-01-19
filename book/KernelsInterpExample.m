global GAUSSQR_PARAMETERS;
GAUSSQR_PARAMETERS.RANDOM_SEED(0);

% Define data locations and evaluation locations
x = pickpoints(0,1,10,'rand');
xeval = pickpoints(0,1,500);

% Choose a function and evaluate it at those locations
yf = @(x) cos(2*pi*x);
y = yf(x);
yeval = yf(xeval);

% Choose a kernel and shape parameter
rbf = @(e,r) (1+(e*r)).*exp(-(e*r));
ep = 4;

% Form the necessary distance matrices
DM = DistanceMatrix(x,x);
DMeval = DistanceMatrix(xeval,x);

% Compute the kernel matrices from the distance matrices
K = rbf(ep,DM);
Keval = rbf(ep,DMeval);

% Evaluate the interpolant at the desired locations
seval = Keval*(K\y);

% Plot the data, interpolant, and true solution together
plot(x,y,'sk','linewidth',6)
hold on
plot(xeval,seval,'b','linewidth',2)
plot(xeval,yeval,'--r','linewidth',3)
hold off
legend('Data','Interpolant','True Function','location','north')
ylim([-1 1])
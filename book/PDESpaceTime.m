% This script is going to consider the space-time kernel
% We want to solve u_t - heat_const*u_xx = 0
%                  u(0,t) = u(1,t) = 0, u(x,0) = x*(1-x)
% u(x,t) is the solution, which we represent as
%       u(x,t) = sum_{i=1}^N c_i K([x,t],[x_i,t_i])
% where [x_i,t_i] are the points in the space-time domain
% The true solution to this problem is ...
fbc  = @(x) zeros(size(x,1),1);
fic  = @(x) x(:,1).*(1-x(:,1));
fint = @(x) zeros(size(x,1),1);
heat_const = .5;

% We call x both the time and space points
Nx = 35;
Nt = 36;
tmax = 1;
xall = pick2Dpoints([0 0],[1 tmax],[Nx,Nt]);

% Come up with some evaluation points
NNx = 50;
NNt = tmax*NNx;
xeval = pick2Dpoints([0 0],[1 tmax],[NNx,NNt]);

% Choose an RBF for our problem
% epvec(1) is the x shape parameter, epvec(2) is the t shape parameter
epvec = [1,1];
rbfM2 = @(r) (1+r).*exp(-r);
rbfM2dt = @(r,dt,ep_t) -ep_t^2*exp(-r).*dt;
rbfM2dxx = @(r,dx,ep_x) ep_x^2*exp(-r).*(ep_x*dx.^2./(r+eps)-1);

% We organize the points so that:
%     [x=0]
%     [x=1]
%     [t=0]
%     [interior]
xbc  = [xall(xall(:,1)==0,:);xall(xall(:,1)==1,:)];
xic  = xall(xall(:,2)==0 & xall(:,1)~=0 & xall(:,1)~=1,:);
xint = xall(xall(:,2)~=0 & xall(:,1)~=0 & xall(:,1)~=1,:);
x = [xbc;xic;xint];

% Evaluate the collocation matrix
Abc  = rbfM2(DistanceMatrix(xbc,x,epvec));
Aic  = rbfM2(DistanceMatrix(xic,x,epvec));
Aint = rbfM2dt(DistanceMatrix(xint,x,epvec),DifferenceMatrix(xint(:,2),x(:,2)),epvec(2)) - ...
       heat_const*rbfM2dxx(DistanceMatrix(xint,x,epvec),DifferenceMatrix(xint(:,1),x(:,1)),epvec(1));
A = [Abc;Aic;Aint];

% Evaluate the RHS
rhsbc  = fbc(xbc);
rhsic  = fic(xic);
rhsint = fint(xint);
rhs = [rhsbc;rhsic;rhsint];

% Solve the system to find the coefficients
coef = A\rhs;

% Evaluate at points in the domain
yeval = rbfM2(DistanceMatrix(xeval,x,epvec))*coef;

% Plot the results
X = reshape(xeval(:,1),NNx,NNt);
T = reshape(xeval(:,2),NNx,NNt);
Y = reshape(yeval,NNx,NNt);
surf(X,T,Y,'edgecolor','none');
xlabel('x')
ylabel('t')
zlabel('u')
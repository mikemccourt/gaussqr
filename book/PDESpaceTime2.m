% This script is going to consider the space-time kernel
% Now we want to work with a transport equation
%          u_t + move_const*u_x = 0
%          u(-1,t) = 0, u(x,0) = f(x)
% The true solution is u(x,t) = f(x-c*t)
%
% We approximate u(x,t) as
%       u(x,t) = sum_{i=1}^N c_i K([x,t],[x_i,t_i])
% where [x_i,t_i] are the points in the space-time domain
move_const = 2;
fbc  = @(x) zeros(size(x,1),1);
fic  = @(x) (64*(-x(:,1)).^3.*(1+(x(:,1))).^3).*(x(:,1)<=0);
true_sol  = @(x) (64*(-(x(:,1)-move_const*x(:,2))).^3.*(1+(x(:,1)-move_const*x(:,2))).^3).*(abs(x(:,1)-move_const*x(:,2)+.5)<=.5);

% We call x both the time and space points
Nx = 35;
Nt = 36;
tmax = 1;
xall = pick2Dpoints([-1 0],[1 tmax],[Nx,Nt]);

% Come up with some evaluation points
NNx = 50;
NNt = tmax*NNx;
xeval = pick2Dpoints([-1 0],[1 tmax],[NNx,NNt]);

% Choose to evaluate error at evaluation points or collocation points
%      1 - evaluation points
%      2 - collocation points
err_eval = 1;

% Choose an RBF for our problem
% epvec(1) is the x shape parameter, epvec(2) is the t shape parameter
% Note rbfM2dt and rbfM2dx are the same, I am just distinguishing for
% clarity because I don't want to make a silly mistake
epvec = [1,1];
rbfM2 = @(r) (1+r).*exp(-r);
rbfM2dt = @(r,dt,ep_t) -ep_t^2*exp(-r).*dt;
rbfM2dx = @(r,dx,ep_x) -ep_x^2*exp(-r).*dx;

% We organize the points so that:
%     [x=0]
%     [x=1]
%     [t=0]
%     [interior]
xbc  = xall(xall(:,1)==-1,:);
xic  = xall(xall(:,2)==0 & xall(:,1)~=-1,:);
xint = xall(xall(:,2)~=0 & xall(:,1)~=-1,:);
x = [xbc;xic;xint];

% Evaluate the collocation matrix
Abc  = rbfM2(DistanceMatrix(xbc,x,epvec));
Aic  = rbfM2(DistanceMatrix(xic,x,epvec));
Aint = rbfM2dt(DistanceMatrix(xint,x,epvec),DifferenceMatrix(xint(:,2),x(:,2)),epvec(2)) + ...
       move_const*rbfM2dx(DistanceMatrix(xint,x,epvec),DifferenceMatrix(xint(:,1),x(:,1)),epvec(1));
A = [Abc;Aic;Aint];

% Evaluate the RHS
rhsbc  = fbc(xbc);
rhsic  = fic(xic);
rhsint = zeros(size(xint,1),1);
rhs = [rhsbc;rhsic;rhsint];

% Solve the system to find the coefficients
coef = A\rhs;

% Evaluate at points in the domain
ueval = rbfM2(DistanceMatrix(xeval,x,epvec))*coef;

% Plot the results
h = figure;
subplot(1,2,1)
X = reshape(xeval(:,1),NNx,NNt);
T = reshape(xeval(:,2),NNx,NNt);
U = reshape(ueval,NNx,NNt);
surf(X,T,U,'edgecolor','none');
xlabel('x')
ylabel('t')
zlabel('u')

% Determine which points to check the error at
if err_eval==1
    xtest = xeval;
    Nxtest = NNx;
    Nttest = NNt;
else
    xtest = xall;
    Nxtest = Nx;
    Nttest = Nt;
end

% Compare this solution to the true solution at all the collocation points
utest = rbfM2(DistanceMatrix(xtest,x,epvec))*coef;
utrue = true_sol(xtest);
title(sprintf('Relative error %g',errcompute(utest,utrue)))

% Plot the difference between the collocation and
% Normalize it against the norm of the initial condition
subplot(1,2,2)
udiff = abs(utest-utrue)/norm(utrue(:,1));
X = reshape(xtest(:,1),Nxtest,Nttest);
T = reshape(xtest(:,2),Nxtest,Nttest);
Udiff = reshape(udiff,Nxtest,Nttest);
surf(X,T,Udiff,'edgecolor','none')
xlabel('x')
ylabel('t')
zlabel('relative error')

% Compare to the best possible error that could occur, which would occur
% for the interpolation without any PDE stuff at all
coef_int = rbfM2(DistanceMatrix(x,x,epvec))\true_sol(x);
uint = rbfM2(DistanceMatrix(xtest,x,epvec))*coef_int;

h_int = figure;
subplot(1,2,1)
Uint = reshape(uint,Nxtest,Nttest);
surf(X,T,Uint,'edgecolor','none')
xlabel('x')
ylabel('t')
zlabel('interpolation solution')
title(sprintf('Relative error %g',errcompute(uint,utrue)))

subplot(1,2,2)
uintdiff = abs(uint-utrue)/norm(utrue(:,1));
Uintdiff = reshape(uintdiff,Nxtest,Nttest);
surf(X,T,Uintdiff,'edgecolor','none')
xlabel('x')
ylabel('t')
zlabel('interpolation error')
% This script is going to consider the space-time kernel
% We want to solve u_t - heat_const*u_xx = 0
%                  u(0,t) = u(1,t) = 0, u(x,0) = x*(1-x)
% u(x,t) is the solution, which we represent as
%       u(x,t) = sum_{i=1}^N c_i K([x,t],[x_i,t_i])
% where [x_i,t_i] are the points in the space-time domain

global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

fbc  = @(x) zeros(size(x,1),1);
fic  = @(x) 4*x(:,1).*(1-x(:,1));
fint = @(x) zeros(size(x,1),1);
heat_const = .25;

% We call x both the time and space points
Nx = 15;  Nt = 20;  tmax = 1;
% Create halton points on the interior, uniform on boundary
xint = pick2Dpoints([0 0],[1 tmax],[Nx-2,Nt-1],'halt');
xbc = [kron([0;1],ones(Nt,1)),repmat(pickpoints(tmax/Nt,tmax,Nt),2,1)];
xic = [pickpoints(0,1,Nx),zeros(Nx,1)];
x = [xbc;xic;xint];

% Come up with some evaluation points
NNx = 50;
NNt = tmax*NNx;
xeval = pick2Dpoints([0 0],[1 tmax],[NNx,NNt]);

% Choose a kernel for our problem
% epvec(1) is the x shape parameter, epvec(2) is the t shape parameter
epvec = [.5,1];
rbfM2 = @(r) (1+r).*exp(-r);
rbfM2dt = @(r,dt,ep_t) -ep_t*exp(-r).*(ep_t*dt);
rbfM2dxx = @(r,dx,ep_x) ep_x^2*exp(-r).*((ep_x*dx).^2./(r+eps)-1);

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

% Solve the system to find the solution coefficients
coef = A\rhs;

% Evaluate at points in the domain
ueval = rbfM2(DistanceMatrix(xeval,x,epvec))*coef;
toc

% Plot the results
X = reshape(xeval(:,1),NNx,NNt);
T = reshape(xeval(:,2),NNx,NNt);
U = reshape(ueval,NNx,NNt);
surf(X,T,U,'edgecolor','none');
xlabel('x')
ylabel('t')
zlabel('u')
view([1 1 1])
zlim([-.002,1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a "true" solution that is based on finite differences
% This must be at a significantly finer resolution to work
tic
FD_res = 100;
Nxfd = FD_res*Nx;
delta_x = 1/(Nxfd-1);
x1d = pickpoints(0,1,Nxfd);
Ntfd = FD_res*Nt;
delta_t = 1/Ntfd;

% Initialize data from problem
uFD = zeros(Nxfd,Ntfd+1);
uFD(:,1) = fic(x1d);
t = 0;

% Create the finite difference second derivative operator
% Zero out the first and last rows since they are BC
FDmat = -gallery('tridiag',Nxfd)/delta_x^2;
FDmat([1,end],:) = zeros(2,size(FDmat,2));

% Choose a time-stepping scheme
%    1 - backward Euler
%    2 - 2nd order BDF
time_scheme = 2;

% Create the backward Euler solution operator
BEmat = speye(Nxfd,Nxfd) - delta_t*heat_const*FDmat;

if time_scheme==1
    % Perform the time stepping with backward Euler
    for k=1:Ntfd
        t = t + delta_t;
        uFD(:,k+1) = BEmat\uFD(:,k);
    end
elseif time_scheme==2
    % Do one backward Euler step to get going
    uFD(:,2) = BEmat\uFD(:,1);
    
    % Create the 2nd order BDF operator
    BDFmat = speye(Nxfd,Nxfd) - 2/3*delta_t*heat_const*FDmat;
    
    % Perform the time stepping with 2nd order BDF
    for k=1:Ntfd-1
        t = t + delta_t;
        xt = [x1d,t*ones(size(x1d,1),1)];
        rhs = 4/3*uFD(:,k+1) - 1/3*uFD(:,k) + 2/3*delta_t*fint(xt);
        uFD(:,k+2) = BDFmat\rhs;
    end
end
% This is for forward Euler, which is unwise
%     uFD(:,k+1) = uFD(:,k) + delta_t*heat_const*FDmat*uFD(:,k);

uFD_final = uFD(:,end);
toc

% Compare the solution from finite differences to the space-time
% collocation solution
xtest = [x1d,tmax*ones(size(x1d,1),1)];
utest = rbfM2(DistanceMatrix(xtest,x,epvec))*coef;
fprintf('Relative difference, %d RBF, %d FD: %g\n',Nx,Nxfd,errcompute(utest,uFD_final))

% Plot the collocation points
h_pts = figure;
plot(xint(:,1),xint(:,2),'.r','markersize',10)
hold on
plot(xbc(:,1),xbc(:,2),'ob','linewidth',2)
plot(xic(:,1),xic(:,2),'+k','linewidth',2)
hold off
xlabel('x')
ylabel('t')
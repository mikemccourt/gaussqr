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

% Solve the system to find the solution coefficients
coef = A\rhs;

% Evaluate at points in the domain
ueval = rbfM2(DistanceMatrix(xeval,x,epvec))*coef;

% Plot the results
X = reshape(xeval(:,1),NNx,NNt);
T = reshape(xeval(:,2),NNx,NNt);
U = reshape(ueval,NNx,NNt);
surf(X,T,U,'edgecolor','none');
xlabel('x')
ylabel('t')
zlabel('u')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a "true" solution that is based on finite differences
% This must be at a significantly finer resolution to work
FD_res = 50;
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

% Compare the solution from finite differences to the space-time
% collocation solution
xtest = [x1d,tmax*ones(size(x1d,1),1)];
utest = rbfM2(DistanceMatrix(xtest,x,epvec))*coef;
fprintf('Relative difference, %d RBF, %d FD: %g\n',Nx,Nxfd,errcompute(utest,uFD_final))
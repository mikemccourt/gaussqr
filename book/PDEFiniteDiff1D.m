% PDEFiniteDiff1D
% Here we solve a first-order variable coefficient wave equation
%     u_t + c(x)u_x = 0, x = [0,2*pi], t = [0,8]
%                        c(x) = .2 + sin(x-1)^2
%                        u(x,0) = exp(-100*(x-1)^2)
% We consider N=400 uniform discretization points and use ode45 to conduct
% the time stepping explicitly.  There aren't really any boundary
% conditions needed for this problem, so we specify none.

% Choose our physical discretization: just a uniform grid although it could
% be scattered if we wanted it to be.  Note that these points are also the
% points on which the solution will be approximated; evaluation at other
% points would require interpolation of the solution.
N = 200;
x = pickpoints(0,2*pi,N);
Nx = 13;

% Choose our time span and max allowed time step
tspan = [0 8];
odeopts = odeset('MaxStep',.0125);

% Define the wave speed coefficients
c = .2 + sin(x-1).^2;

% Define the initial condition
u0 = exp(-100*(x-1).^2);

% Choose the C4 Matern kernel to approximate spatial derivatives
rbf = @(e,r) (1+e*r+(e*r).^2/3).*exp(-e*r);
rbfx = @(e,r,dx) -e^2/3*dx.*(1+e*r).*exp(-e*r);
ep = 1;

%%%%%%
% Create the finite difference operator
% Find the nearest neighbors to the solution points
nearest = num2cell(knnsearch(x,x,'K',Nx),2);
% Find the finite difference coefficients
FDcell = cellfun(@(xe,xi) ...
    rbfx(ep,DistanceMatrix(xe,x(xi,:)),DifferenceMatrix(xe,x(xi,:)))/...
    rbf(ep,DistanceMatrix(x(xi,:),x(xi,:))),...
    num2cell(x,2),nearest,'UniformOutput',0);
% Form the vectors with the values for the sparse matrix
FDvecs = cell2mat(cellfun(@(row,cols,vals)[row*ones(1,Nx);cols;vals],...
    num2cell((1:N)',2),nearest,FDcell,'UniformOutput',0)');
% Create the sparse matrix from those vectors
LFD = sparse(FDvecs(1,:),FDvecs(2,:),FDvecs(3,:),N,N);
%%%%%%

% Change the first row to a boundary condition
LFD(1,:) = zeros(1,N);

% Define a function that acts like u_t = f(t,u)
% In this case, f(t,u) = -c(x)*u_x
odefun = @(t,u) -bsxfun(@times,c,LFD*u);

% Call the ODE solver
[T,U] = ode15s(odefun,tspan,u0,odeopts);

h = figure;
surf(x,T,U,'edgecolor','none'),view([0 0 1]),grid off
xlim([0,2*pi]),xlabel('x'),ylabel('t')
colormap gray,C = colormap;colormap(flipud(C))
h_color = colorbar('limits',[-5.5e-04,1],'ticks',[0 .5 1]);
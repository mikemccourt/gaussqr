% ex12b.m
% This example solves the Burgers' equation
%   u_t + u*u_x = u_xx/R
% with u(0,t)=u(1,t)=0, and u(x,0) = sin(pi*x)
% The purpose of this is to compare to Hon's paper
%
% We're using their semi-implicit time discretization
%  u^m + dt(u^{m-1}u_x^m - 1/R u_{xx}^m) = u^{m-1}
%
% Because N=11 for this problem, we'll use GaussQR
%
% NOTE: In a first attempt to make this work, I have
% changed the domain of the problem to [-1,1].  We'll
% see if this is necessary - maybe switch to [-1/2,1/2]
% Because of this, the errors cannot be computed right now
% See: An efficient numerical scheme for Burgers' equation
% YC Hon, XZ Mao - Applied Mathematics and Computation, 1998
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION = 1;

N = 13;
dt = .001;
T = 1;
R = 100;
t_steps = dt:dt:T;
coefs = zeros(N,length(t_steps)+1);
I = eye(N);
output_frequency = 2;

x = pickpoints(-1,1,N);
IC = @(x) sin(pi*x/2+pi/2);
uold = IC(x);

alpha = 1;

% Do an epsilon search based on IC
epsearch = 1;

if epsearch
    NN = 100;
    xx = pickpoints(-1,1,NN);
    uu = IC(xx);
    ep = fminbnd(@(ep)errcompute(gqr_eval(gqr_solve(x,uold,ep,alpha),xx),uu),.1,1);
else
    ep = .1;
end

% Here we interpolate the initial condition as evaluated at the collocation points
% This will allow us to create an initial set of coefficients for use as the initial condition
% errs = [];
GQR = gqr_solve(x,uold,ep,alpha); % Set up GQR object
u = gqr_eval(GQR,x);
coefs(:,1) = GQR.coef;
% errs(1) = errcompute(u,uold);

% Set up the Phi basis which we use for the collocation within the Psi basis
% Notice that we can reuse the Phi values from earlier, as well as the Rbar values
% The derivatives must be computed here
Phi = [GQR.stored_phi1,GQR.stored_phi2];
Rbar = GQR.Rbar;
Phix = gqr_phi(GQR,x,1);
Phixx = gqr_phi(GQR,x,2);

% We can form the linear portion of the solution matrix here, which will be added
% to the nonlinear portion (dependent on last time step data) at each time step
% Notice also that this includes the boundary terms which will be overwritten when
% forming the full collocation solution matrix
Lp = 1/dt*Phi - 1/R*Phixx;

h = figure;

k = 2; % IC takes up first storage spot
for t=t_steps
	% Add the nonlinear terms to the matrix of linear terms
	% Note that the nonlinear terms are a function of the previous solution
    A = Lp + diag(uold)*Phix;
	
	% Replace the top and bottom rows with the identity operator for the BC
	% The identity is the basis with no derivatives in it so that u(-1,t) = u(1,t) = 0
    A([1,end],:) = Phi([1,end],:); % Impose BC
    
    % Form the collocation matrix in the stable basis
    A = A*[I;Rbar];
    
    % Create the right hand side, using both collocation and BC
    rhs = 1/dt*uold;
    rhs([1,end]) = [0;0];
    
    % Solve the system to determine the solution coefficients
    % Evaluate the solution and store the coefficients
    GQR.coef = A\rhs;
    u = gqr_eval(GQR,x);
    coefs(:,k) = GQR.coef;
    
    % Compute the error of the solution ... not working yet
    
    uold = u;
    k = k+1;
    if mod(k,output_frequency)==0 
        plot(x,u,'linewidth',3)
        ylim([0 1])
        title(sprintf('time %4.3f',t))
        pause(.00001)
    end
end
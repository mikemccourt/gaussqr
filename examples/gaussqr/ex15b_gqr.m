% ex15b.m
% This example solves the Burgers' equation
%   u_t + u*u_x = u_xx/R
% with u(0,t)=u(1,t)=0, and u(x,0) = sin(pi*x)
% The purpose of this is to compare to Hon's paper
%
% We're using their semi-implicit time discretization
%  u^m + dt(u^{m-1}u_x^m - 1/R u_{xx}^m) = u^{m-1}
%
% Because N=11 for this problem, we'll use GaussQR
rbfsetup
global GAUSSQR_PARAMETERS

N = 13;
dt = .001;
T = 1;
R = 100;
t_steps = dt:dt:T;
coefs = zeros(N,length(t_steps)+1);
I = eye(N);

x = pickpoints(0,1,N);
IC = @(x) sin(pi*x);
uold = IC(x);

alpha = 1;

% Do an epsilon search based on IC
epsearch = 0;
% Output time along the way
output_time = 0;

if epsearch
    NN = 100;
    xx = pickpoints(0,1,NN);
    uu = IC(xx);
    epvec = logspace(-1,1,55);
    errvec = [];
    k = 1;
    for ep=epvec
        GQR = gqr_solve(x,uold,ep,alpha);
        errvec(k) = errcompute(gqr_eval(GQR,xx),uu);
        k = k+1;
    end
    ep = fminbnd(@(ep)errcompute(gqr_eval(gqr_solve(x,uold,ep,alpha),xx),uu),.1,1);
else
    ep = .1;
end

% Need to set up solver
errs = [];
GQR = gqr_solve(x,uold,ep,alpha); % Set up GQR object
[u,GQR] = gqr_eval(GQR,x); % Store GQR x values for later
coefs(:,1) = GQR.coef;
errs(1) = errcompute(u,uold);

% Set up the psi basis
[ep,alpha,Marr,lam] = gqr_solveprep(0,x,ep,alpha);
phi = gqr_phi(Marr,x,ep,alpha);
phi1 = gqr_phi(Marr,x,ep,alpha,1);
phi2 = gqr_phi(Marr,x,ep,alpha,2);
Rbar = GQR.Rbar;
Lp = phi-dt/R*phi2;

k = 2; % IC takes up first storage spot
for t=t_steps
    A = (Lp + dt*diag(uold)*phi1); % Impose PDE
    A([1,end],:) = phi([1,end],:); % Impose BC
    A = A*[I;Rbar]; % RBF-QR
    rhs = uold;
    rhs([1,end]) = [0;0];
    
    GQR.coef = A\rhs;
    u = gqr_eval(GQR,x);
    coefs(:,k) = GQR.coef;
    errs(k) = errcompute(u,ex15b_gqr_truesol(x,t,R));
    
    uold = u;
    k = k+1;
    if mod(k,10)==0 & output_time
        fprintf('%g\n',t)
    end
end
% StableCollDemo.m
% Solves a 2-pt BVP using the stable basis
% Calls on: gqr_solveprep, gqr_phi, gqr_eval, errcompute, pickpoints
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = 100;
ep_vals = logspace(-1,1,30); alpha = .75;
N = 25; x = pickpoints(0,pi,N);
N_eval = 100; x_eval = pickpoints(0,pi,N_eval);
f = @(x) -sin(x); sol = @(x) sin(x);
u_eval = sol(x_eval);
err_vals = zeros(size(ep_vals));
k = 1;
for ep=ep_vals
    GQR = gqr_solveprep(0,x,ep,alpha); Rbar = GQR.Rbar;
    Phi_interior = gqr_phi(GQR,x(2:end-1),2);
    Phi_boundary = gqr_phi(GQR,x([1,end]));
    Psi = [Phi_interior;Phi_boundary]*[eye(N);Rbar];
    rhs_interior = f(x(2:end-1));
    rhs_boundary = [0;0];
    rhs = [rhs_interior;rhs_boundary];
    GQR.coef = Psi\rhs;
    u_gaussqr = gqr_eval(GQR,x_eval);
    err_vals(k) = errcompute(u_gaussqr,u_eval);
    k = k + 1;
end
loglog(ep_vals,err_vals,'linewidth',3)
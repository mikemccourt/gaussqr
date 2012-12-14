% UnstableBVPDemo.m
% This shows the ill-conditioning of standard kernel collocation
% Calls on : DistanceMatrix, pickpoints
rbf = @(e,r) exp(-(e*r).^2);
d2rbf = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);
ep = .1;
f = @(x) -sin(x); sol = @(x) sin(x);
x = pickpoints(0,pi,25);
DM_interior = DistanceMatrix(x(2:end-1),x);
A_interior = d2rbf(ep,DM_interior);
DM_boundary = DistanceMatrix(x([1,end]),x);
A_boundary = rbf(ep,DM_boundary);
rhs_interior = f(x(2:end-1));
rhs_boundary = [0;0];
A = [A_interior;A_boundary];
rhs = [rhs_interior;rhs_boundary];
coef = A\rhs;
x_eval = pickpoints(0,pi,200);
DM_eval = DistanceMatrix(x_eval,x);
A_eval = rbf(ep,DM_eval);
u_eval = A_eval*coef;
plot(x_eval,u_eval,'b','linewidth',2),xlim([0,pi]),ylim([0,1])
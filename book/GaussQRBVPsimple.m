% GaussQRBVPsimple.m
% This solves a heat equation using GaussQR
% The steady state solution is u(x) = .5*x*(1-x)
% Calls on: pickpoints, gqr_solveprep, gqr_phi, gqr_eval
N = 15;N_eval = 100;dt = .01;T = 1;ep = 1e-4;alpha = 2;
tvals = (0:dt:T);sols = zeros(N_eval,length(tvals));
% Set up domain and GQR object
x = pickpoints(0,1,N);
GQR = gqr_solveprep(0,x,ep,alpha);
x_eval = pickpoints(0,1,100);
% Identify boundary and interior (nonboundary) points
xBC = x(find((x==0) + (x==1))); xINT = setdiff(x,xBC);
% Set up needed matrices
IRbar = [eye(N);GQR.Rbar];
PhiINT = gqr_phi(GQR,xINT);
allzeros = zeros(length(xBC),length(GQR.Marr));
PhiL = gqr_phi(GQR,xINT,2);
PhiB = gqr_phi(GQR,xBC);
Psi0 = [PhiINT;allzeros]*IRbar;
PsiLB = [PhiL;PhiB]*IRbar;
A = Psi0 - dt*PsiLB;
f = ones(length(xINT),1);
g = zeros(length(xBC),1);
h = [f;g];
% Find initial conditions
PhiIC = gqr_phi(GQR,x);initvals = zeros(N,1);
b = (PhiIC*IRbar)\initvals; GQR.coef = b;
sols(:,1) = gqr_eval(GQR,x_eval);
% Time step to completion
k = 2;
for t = tvals(2:end)
    b = A\(Psi0*b + dt*h); GQR.coef = b;
    sols(:,k) = gqr_eval(GQR,x_eval);
    k = k + 1;
end
[XX,TT] = meshgrid(x_eval,tvals);surf(XX,TT,sols')
shading interp,colormap gray,xlabel('x'),ylabel('t')
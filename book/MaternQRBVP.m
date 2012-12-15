% MaternQRBVP.m
% This solves a heat equation using MaternQR
% Calls on: pickpoints, mqr_solveprep, mqr_phi, mqr_eval
N = 15;dt = .01;T = 1;ep = 1;beta = 6;L = 1;
tvals = (0:dt:T);sols = zeros(N_eval,length(tvals));
% Set up domain and MQR object
x = pickpoints(0,L,N+2);
x = x(2:end-1); %ignore boundary points
MQR = mqr_solveprep(x,L,ep,beta);
x_eval = pickpoints(0,L,100);
% Set up needed matrices
IRbar = [eye(N);MQR.Rbar];
PhiINT = mqr_phi(MQR,xINT);
allzeros = zeros(length(xBC),length(MQR.Marr));
PhiL = mqr_phi(MQR,xINT,2);
PhiB = mqr_phi(MQR,xBC);
Psi0 = [PhiINT;allzeros]*IRbar;
PsiLB = [PhiL;PhiB]*IRbar;
A = Psi0 - dt*PsiLB;
f = ones(length(xINT),1);
g = zeros(length(xBC),1);
h = [f;g];
% Find initial conditions
PhiIC = mqr_phi(MQR,x);initvals = zeros(N,1);
b = (PhiIC*IRbar)\initvals; MQR.coef = b;
sols(:,1) = mqr_eval(MQR,x_eval);
% Time step to completion
k = 2;
for t = tvals(2:end)
    b = A\(Psi0*b + dt*h); MQR.coef = b;
    sols(:,k) = mqr_eval(MQR,x_eval);
    k = k + 1;
end
[XX,TT] = meshgrid(x_eval,tvals);surf(XX,TT,sols')
shading interp,colormap gray,xlabel('x'),ylabel('t')
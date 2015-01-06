% PlotMaternKernelSeries
% Computes and plots kernels given by a Mercer series
close all
x = linspace(0,1,11)';
xx = linspace(0,4,1201)';
%% C0 Matern kernel on [0,\infty), special case \alpha = 4\epsilon
ep = 1;
phifunc = @(n,x) sin(exp(-2*ep*x)*(2*n-1)*pi/2)./(exp(-ep*x)*ones(size(n)));
lambdafunc = @(n) 8./((2*n-1)*pi).^2;
M = 100;
%% Mercer series
N = length(x);
Lambda = diag(lambdafunc(1:M));
Phi_interp = phifunc(1:M,x);
Phi_eval = phifunc(1:M,xx);
Kbasis = Phi_eval*Lambda*Phi_interp';
%% Plot kernel basis obtained via Mercer series
plot(xx,Kbasis)
%% Plot first 4 eigenfunctions
figure
plot(xx,Phi_eval(:,1:4))
%% HS-SVD stuff
I_N = eye(N);
Phi_1 = Phi_interp(:,1:N);
Phi_2 = Phi_interp(:,N+1:end);
Lambda_1 = Lambda(1:N,1:N);
Lambda_2 = Lambda(N+1:M,N+1:M);
Correction = Lambda_2*(Phi_1\Phi_2)'/Lambda_1;
Psi_basis = Phi_eval*[I_N;Correction];
%Kbasis = Psi_eval*Lambda_1*Phi_1'/Lambda(1,1);
%% Plot Psi basis
figure
plot(xx,Psi_basis);%(:,1:3))
%% Plot correction added to first N eigenfunctions to obtain Psi basis
Phi_correct = Phi_eval(:,N+1:M)*Correction;
figure
plot(xx,Phi_correct(:,1:5))


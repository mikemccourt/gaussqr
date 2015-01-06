% PlotKernelSeries
% Computes and plots kernels given by a Mercer series
close all
x = linspace(0,1,11)';
xx = linspace(0,1,1201)';
%% iterated Brownian bridge kernel
x=x(2:end-1);
ep = 50; beta = 1;
phifunc = @(n,x) sqrt(2)*sin(pi*x*n);
lambdafunc = @(n) ((n*pi).^2+ep^2).^(-beta);
if beta < 3
    M = 1000;
else
    M = ceil(1/pi*sqrt(eps^(-1/beta)*(N^2*pi^2+ep^2)-ep^2));
end
%% Mercer series
N = length(x);
Lambda = diag(lambdafunc(1:M));
Phi_interp = phifunc(1:M,x);
Phi_eval = phifunc(1:M,xx);
Kbasis = Phi_eval*Lambda*Phi_interp'/Lambda(1,1);
%% Plot cardinal basis obtained via Mercer series
plot(xx,Kbasis)
%% Plot first 5 eigenfunctions
figure
plot(xx,Phi_eval(:,1:5))
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
Phi_correct = Phi_eval(:,N+1:M)*Correction;
%% Plot correction added to first N eigenfunctions to obtain Psi basis
figure
plot(xx,Phi_correct(:,1:5))


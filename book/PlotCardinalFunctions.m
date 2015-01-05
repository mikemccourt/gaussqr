% PlotKernelSeries
% Computes and plots kernels given by a Mercer series
close all
x = linspace(0,1,22)';
%x = pickpoints(0,1,22,'cheb');
xx = linspace(0,1,1201)';
%% iterated Brownian bridge kernel
x=x(2:end-1);
ep = 50; beta = 20;
phifunc = @(n,x) sqrt(2)*sin(pi*x*n);
lambdafunc = @(n) ((n*pi).^2+ep^2).^(-beta);
if beta < 3
    M = 1000;
else
    M = ceil(1/pi*sqrt(eps^(-1/beta)*(N^2*pi^2+ep^2)-ep^2))
end
%% Mercer series
N = length(x);
Lambda = diag(lambdafunc(1:M));
Phi_interp = phifunc(1:M,x);
Phi_eval = phifunc(1:M,xx);
Kmatrix = Phi_interp*Lambda*Phi_interp';
Kbasis = Phi_eval*Lambda*Phi_interp';
cardinalbasis = (Kmatrix\Kbasis')';
%% Plot cardinal basis obtained via Mercer series
plot(xx,cardinalbasis(:,[1 5 10 15 20]))
%% HS-SVD stuff
I_N = eye(N);
Phi_1 = Phi_interp(:,1:N);
Phi_2 = Phi_interp(:,N+1:end);
Lambda_1 = Lambda(1:N,1:N);
Lambda_2 = Lambda(N+1:M,N+1:M);
Correction = Lambda_2*(Phi_1\Phi_2)'/Lambda_1;
Psi_matrix = Phi_interp*[I_N;Correction];
Psi_basis = Phi_eval*[I_N;Correction];
cardinalbasis = (Psi_matrix\Psi_basis')';
%% Plot cardinal basis obtained via HS-SVD
figure
plot(xx,cardinalbasis(:,[1 5 10 15 20]))


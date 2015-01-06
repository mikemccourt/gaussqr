% IBBKernelEx1.m
% Computes and plots iterated Brownian bridge kernels via their Mercer series
x = linspace(0,1,11)'; xx = linspace(0,1,1201)';
%% iterated Brownian bridge kernel
x=x(2:end-1); N = length(x);
ep = 50; beta = 20;
phifunc = @(n,x) sqrt(2)*sin(pi*x*n);
lambdafunc = @(n) ((n*pi).^2+ep^2).^(-beta);
%% Mercer series
if beta < 3 
    M = 1000;
else
    M = ceil(1/pi*sqrt(eps^(-1/beta)*(N^2*pi^2+ep^2)-ep^2));
end 
Lambda = diag(lambdafunc(1:M));
Phi_interp = phifunc(1:M,x);
Phi_eval = phifunc(1:M,xx);
Kbasis = Phi_eval*Lambda*Phi_interp'/Lambda(1,1);
%% Plot kernel basis obtained via Mercer series
plot(xx,Kbasis)

% IBBKernelEx2.m
% Computes and plots cardinal functions for 
% iterated Brownian bridge and Gaussian kernels
N = 20; NN = 300; pts = {'even','cheb'};
xeval = pickpoints(0,1,NN);
%% Iterated Brownian bridge kernel
ep = 50; beta = 20; 
phifunc = @(n,x) sqrt(2)*sin(pi*x*n);       % eigenfunctions
% Truncation length
M = ceil(1/pi*sqrt(eps^(-1/beta)*(pi^2+ep^2)-ep^2));
n = 1:M;
lambda = ((pi*n).^2+ep^2).^(-beta);         % eigenvalues
neval = [1 5 10 15 20]; % index of cardinal functions for plot
for k=1:2
    x = pickpoints(0,1,N+2,pts{k}); x = x(2:end-1);
    Phi = phifunc(n,x);             % eigenfunction matrix for interpolation 
    Phieval = phifunc(n,xeval);           %            and for evaluation
    % Hilbert-Schmidt SVD
    [Q,R] = qr(Phi);
    R1 = R(:,1:N);    R2 = R(:,N+1:end);
    opts.UT = true;   Rhat = linsolve(R1,R2,opts);
    D1 = repmat(lambda(1:N),M-N,1);
    D2 = repmat(lambda(N+1:end)',1,N);
    Correction = D2.*Rhat'./D1;
    Psi = Phi*[eye(N);Correction];          % stable interpolation matrix
    Psieval = Phieval*[eye(N);Correction];  %  and for evaluation
    % Evaluate the cardinal functions and plot them
    cardfuncs = Psieval/Psi; 
    subplot(2,2,k), plot(xeval,cardfuncs(:,neval),'linewidth',3)
end  
%% Gaussian kernel
rbf = @(e,r) exp(-(e*r).^2);  
ep = 5.75;                                  % to match scaling of K_{20,50}
neval = neval+1;          % small shift since Gaussians nonzero at boundary
for k=1:2
    x = pickpoints(0,1,N+2,pts{k});
    DM = DistanceMatrix(x,x);           % distance matrix for interpolation
    DMeval = DistanceMatrix(xeval,x);         %       and for evaluation
    K = rbf(ep,DM);                     % interpolation matrix
    Keval = rbf(ep,DMeval);             %   and for evaluation
    % Evaluate cardinal functions and plot them 
    cardfuncs = Keval/K; 
    subplot(2,2,k+2), plot(xeval,cardfuncs(:,neval),'linewidth',3)
end 

% PlotCardinalFunctions
% This code creates figure 3 from the paper with Roberto
% The pictures in it were nontrivial to create, so the code below may look
% kind of odd or poorly written

% Choose evenly spaced points to compute with
% Remember that we dump the endpoints for the Matern interpolation
N = 20; NN = 300;
xeval = pickpoints(0,1,NN);
pts = {'even','cheb'}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Iterated Brownian bridge first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choice of kernel as given in the paper
ep = 50; beta = 20;

% Eigenfunctions
% The sqrt(2) is unnecessary but is present for normalization purposes
phifunc = @(n,x) sqrt(2)*sin(pi*x*n);

% Truncation length
M = ceil(1/pi*sqrt(eps^(-1/beta)*(pi^2+ep^2)-ep^2));
n = 1:M;

% Eigenvalues
lambda = ((pi*n).^2+ep^2).^(-beta);

% Choose which cardinal functions we want to evaluate
neval = [1 5 10 15 20];

for k=1:2
    % Choose the points for computing
    x = pickpoints(0,1,N+2,pts{k});
    x = x(2:end-1);
    
    % Fill the eigenfunction matrix for interpolation and evaluation
    Phi = phifunc(n,x);
    Phieval = phifunc(n,xeval);

    % Compute the QR decomposition of the short, fat matrix Phi
    % This is also probably unneeded, but whatever
    [Q,R] = qr(Phi);
    R1 = R(:,1:N);    R2 = R(:,N+1:end);

    % Apply inv(R1)*R2
    opts.UT = true;
    Rhat = linsolve(R1,R2,opts);

    % Compute the eigenvalue matrix
    D1 = repmat(lambda(1:N),M-N,1);
    D2 = repmat(lambda(N+1:end)',1,N);

    % Form the Rbar correction matrix
    Correction = D2.*Rhat'./D1;

    % Form the stable basis interpolation matrix
    Psi = Phi*[eye(N);Correction];
    Psieval = Phieval*[eye(N);Correction];

    % Evaluate the cardinal functions and plot them
    cardfuncs = Psieval/Psi;
    subplot(2,2,k)
    plot(xeval,cardfuncs(:,neval),'linewidth',3)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Now Gaussians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the Gaussian
rbf = @(e,r) exp(-(e*r).^2);
ep = 5.75; % As given in the paper

% Small shift to account for non-zeroness of the Gaussian kernels at the
% boundary points: we compute on N=22 points but still want the same
% looking cardinal functions
neval = neval + 1;

for k=1:2
    % Choose the interpolation points and form the distance matrices
    x = pickpoints(0,1,N+2,pts{k});
    DM = DistanceMatrix(x,x);
    DMeval = DistanceMatrix(xeval,x);
    
    % Compute the kernel matrices
    K = rbf(ep,DM);
    Keval = rbf(ep,DMeval);
    
    % Evaluate and plot the cardinal functions
    cardfuncs = Keval/K;
    subplot(2,2,k+2)
    plot(xeval,cardfuncs(:,neval),'linewidth',3)
end
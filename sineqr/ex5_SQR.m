% ex5_SQR
% This puts together convergence order pictures for SineQR
% We consider two choices:
%     General functions with arbitrary "BC" embedded in [0,L]
%     Functions which are homogeneous on [0,L]
% These tests will fix L, sigma and beta and consider an increase in N
% Points will be evenly spaced at first, althought we could consider other
%   choices of point distributions
% For the embedded setting, the embedding cushion will be fixed.  The
%   effect of the embedding cushion still needs to be considered.
% I guess we will need to consider both the kernel version and the RBF-QR
%   version of the solve.  Need to think about that because one may be
%   better/faster than the other in different conditions
rbfsetup

% The range of N values to consider
Nvec = 10:5:100;
% The orders (smoothness) of the kernel to consider
betavec = 1:5;
% The kernel shape parameter
sigma = .1;
% The length of the domain
L = 1;
% The embedding width for nonhomogeneous functions
embed_cushion = .1;
% The number of evenly spaced points at which to sample error
NN = 397;

% This determines how many extra basis functions should be added to the
% RBF-QR evaluation to get the necessary accuracy: M = Mfactor*N
% Actually picking a good value for this may be difficult
% I guess the minimum should be something like 1.1
Mfactor = 4.5;

% Define the eigenfunctions and eigenvalues
sinfunc = @(n,L,x) sqrt(2/L)*sin(pi*x*n/L);
lamfunc = @(n,L,sigma,beta) ((pi*n/L).^2+sigma^2).^(-beta);

% This is the function we are interested in considering
% Depending on which function consider, it will choose embedding
fopt = 2;
switch fopt
    case 1
        yf = @(x) sin(2*pi*x) + 1;
        fstr = 'u(x) = sin(2\pi{x})+1, ';
        embed = embed_cushion;
    case 2
        yf = @(x) sin(2*pi*x);
        fstr = 'u(x) = sin(2\pi{x}), ';
        embed = 0;
    otherwise
        error('This function does not exist')
end

% At some point I'll move the solves into a separate function and consider
% warnings there.  For now I'm tired of seeing the warnings pop up.
warning off

errvec = zeros(length(betavec),length(Nvec));
k = 1;
j = 1;
for N=Nvec
    if embed==0
    % Don't consider the homogeneous end points because they are fixed
        [x,spacestr] = pickpoints(0,L,N+2);
        x = x(2:end-1);
    else
        [x,spacestr] = pickpoints(embed*L,(1-embed)*L,N);
    end
    y = yf(x);
    xx = pickpoints(embed*L,(1-embed)*L,NN);
    yy = yf(xx);
    I = eye(N);
    
    j = 1;
    for beta=betavec
        if beta==1 % Work with the kernel form
            K_solve = zeros(N);
            K_eval = zeros(NN,N);
            for n=1:N
                K_solve(:,n) = sobfunc(x,x(n),L,sigma,beta);
                K_eval(:,n) = sobfunc(xx,x(n),L,sigma,beta);
            end
            b = K_solve\y;
            yp = K_eval*b;
        else % Work with the series form
            M = ceil(Mfactor*N);
            n = 1:M;
            S = sinfunc(n,L,x);
            [Q,R] = qr(S);
            R1 = R(:,1:N);
            R2 = R(:,N+1:end);
            opts.UT = true;
            Rhat = linsolve(R1,R2,opts);
            lambda = lamfunc(n,L,sigma,beta);
            D = diag(lambda);
            D1 = diag(lambda(1:N));
            D2 = diag(lambda(N+1:end));
            Rbar = D2*Rhat'/D1;
            b = (S*[I;Rbar])\y;
            SS = sinfunc(n,L,xx);
            yp = (SS*[I;Rbar])*b;
        end

        errvec(j,k) = errcompute(yp,yy);
        j = j + 1;
    end
    
    k = k + 1;
end

warning on

semilogy(Nvec,errvec,'linewidth',2)
xlabel('input points N')
ylabel('RMS relative error')
title(strcat(fstr,sprintf('x\\in[%g,%g], ',embed*L,(1-embed)*L),sprintf('\\sigma=%g',sigma)))
legend('\beta=1','\beta=2','\beta=3','\beta=4','\beta=5','location','northeast')
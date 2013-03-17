% Tests RBF-QR to show that it can avoid the ill-conditioning
% associated with a small shape parameter
% This is considering a 2-pt BVP:
%   u''(x) = f(x)
%   u(0) = u(L) = 0
% I'll need to convert this to look more like the GaussQR stuff
% That'll take time though, and I don't have it right now
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

Nvec = [10,20,40,80,160,320];
betavec = [2,3,4];
NN = 100;
ep = 1;

spaceopt = 'even';

L = 1;
fstr = 'y(x)=exp(x)+x(1-e)-1';
yf = @(x) exp(x)+x*(1-exp(1))-1; % True solution
ff = @(x) exp(x); % Source

yf = @(x) exp(x) - 1 + (1/6).*(x.^3 - 3.*x.^2 + 8.*x - exp(1).*(x.^3+5.*x));
ff = @(x) exp(x) + x*(1-exp(1)) - 1;

% Define the eigenvalues
lamfunc = mqr_solveprep();

% Selecting some points in the domain for computing the error
xx = pickpoints(0,L,NN);
yy = yf(xx);

% Set up the error vectors to store the results
errvec = zeros(length(betavec),length(Nvec));

i = 1;
for beta=betavec
    k = 1;
    for N=Nvec
        K_solve = zeros(N);
        K_eval = zeros(NN,N);
        [x,spacestr] = pickpoints(0,L,N+2,spaceopt);
        x = x(2:end-1); % Dump 0 and L
        y = ff(x);
        I = eye(N);
        
        % Solve it with the QR method
        M = ceil(8.5*N);
        n = 1:M;
        S = mqr_phi(n,x,L);
        [Q,R] = qr(S);
        R1 = R(:,1:N);
        R2 = R(:,N+1:end);
        opts.UT = true;
        Rhat = linsolve(R1,R2,opts);
        lambda = lamfunc(n,L,ep,beta);
        D = diag(lambda);
        D1 = diag(lambda(1:N));
        D2 = diag(lambda(N+1:end));
        Rbar = D2*Rhat'/D1;
        
        S2d = mqr_phi(n,x,L,2);
        M_solve = S2d*[I;Rbar];
        M_eval = mqr_phi(n,xx,L);
        b = M_solve\y;
        yp = M_eval*([I;Rbar]*b);
        errvec(i,k) = errcompute(yp,yy);
        
        k = k+1;
    end
    i = i+1;
end

loglog(Nvec,errvec,'linewidth',2);
xlabel('collocation points')
ylabel('absolute error')

legendvals = [];
for beta=betavec
    legendvals = [legendvals;sprintf('beta=%d',beta)];
end
legend(cellstr(legendvals),'location','southwest')

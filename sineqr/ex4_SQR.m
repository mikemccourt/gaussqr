% Tests RBF-QR to show that it can avoid the ill-conditioning
% associated with a small shape parameter
% This is considering a 2-pt BVP:
%   u''(x) = f(x)
%   u(.1L) and u(.9L) defined
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 4; % Relative RMS

sigmavecd = logspace(0,2,20);
sigmavec =  logspace(0,2,20);
Nvec = [10,20,40];
NN = 100;

spaceopt = 'cheb';
fopt = 'sin';

L = 1;
[yf,fstr] = pickfunc(fopt,1);
fstr = 'y(x)=exp(-x)';
yf = @(x) exp(-x); % True solution
ff = @(x) exp(-x); % Source

% May consider these functions which automatically satisfy the boundary
% conditions f(0)=f(L)=0
% yf = @(x) x.*(L-x);
% yf = @(x) x.*(L-x).*sqrt(x);

% Define the eigenfunctions and eigenvalues
sinfunc = @(n,L,x) sqrt(2/L)*sin(pi*x*n/L);
sinfunc1 = @(n,L,x) (sqrt(2/L)*ones(size(x))*(pi*n/L)).*cos(pi*x*n/L);
sinfunc2 = @(n,L,x) (-sqrt(2/L)*ones(size(x))*(pi*n/L).^2).*sin(pi*x*n/L);
lamfunc = @(n,L,sigma,beta) ((pi*n/L).^2+sigma^2).^(-beta);

% Selecting some points in the domain, not including the boundary
aa = .1*L; bb = .9*L;
xx = pickpoints(aa,bb,NN);
yy = yf(xx);

% Set up the error vectors to store the results
errvec = zeros(length(Nvec),length(sigmavec));
errvecd = zeros(length(Nvec),length(sigmavecd));

% Choose the order of the basis functions
beta = 8;

i = 1;
for N=Nvec
    K_solve = zeros(N);
    K2d = zeros(N);
    K = zeros(N);
    K_eval = zeros(NN,N);
    [x,spacestr] = pickpoints(aa,bb,N,spaceopt);
    y = [yf(x(1));ff(x(2:end-1));yf(x(end))];
    I = eye(N);
    
    k = 1;
    for sigma=sigmavec
        M = ceil(8.5*N);
        n = 1:M;
        S = sinfunc(n,L,x);
        S2d = sinfunc2(n,L,x);
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
        
        Dmat = [S(1,:);S2d(2:end-1,:);S(end,:)]*[I;Rbar];
        b = Dmat\y;
        SS = sinfunc(n,L,xx);
        yp = (SS*[I;Rbar])*b;
        errvec(i,k) = errcompute(yp,yy);
        
        for j=1:N
            K(:,j) = sobfunc(x,x(j),L,sigma,beta);
            K2d(:,j) = sobfunc(x,x(j),L,sigma,beta,2);
            K_eval(:,j) = sobfunc(xx,x(j),L,sigma,beta);
        end
        K_solve = [K(1,:);K2d(2:end-1,:);K(end,:)];
        b = K_solve\y;
        yp = K_eval*b;
        errvecd(i,k) = errcompute(yp,yy);
        k = k+1;
    end
    i = i+1;
end

loglog(sigmavecd,errvecd(1,:),'-bx')
hold on
loglog(sigmavecd,errvecd(2,:),'-g+')
loglog(sigmavecd,errvecd(3,:),'-r^')
loglog(sigmavec,errvec(1,:),'b','LineWidth',3)
loglog(sigmavec,errvec(2,:),'g','LineWidth',3)
loglog(sigmavec,errvec(3,:),'r','LineWidth',3)
hold off
xlabel('\sigma')
ylabel('average error')
ptsstr=strcat(', x\in[',num2str(aa),',',num2str(bb),'],');
title(strcat(fstr,ptsstr,spacestr))
legend('N=10 (Direct)','N=20 (Direct)','N=40 (Direct)','N=10 (QR)','N=20 (QR)','N=40 (QR)', 'Location', 'SouthEast');

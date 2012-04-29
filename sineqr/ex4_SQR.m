% Tests RBF-QR to show that it can avoid the ill-conditioning
% associated with a small shape parameter
% This is considering a 2-pt BVP:
%   u''(x) = f(x)
%   u(.1L) and u(.9L) defined
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 4; % Relative RMS

epvec =  logspace(0,2,20);
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
lamfunc = @(n,L,ep,beta) ((pi*n/L).^2+ep^2).^(-beta);

% Selecting some points in the domain, not including the boundary
aa = .1*L; bb = .9*L;
xx = pickpoints(aa,bb,NN);
yy = yf(xx);

% Set up the error vectors to store the results
errvec = zeros(length(Nvec),length(epvec));
errvecd = zeros(length(Nvec),length(epvec));

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
    for ep=epvec
        M = ceil(8.5*N);
        n = 1:M;
        S = sinfunc(n,L,x);
        S2d = sinfunc2(n,L,x);
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
        
        Dmat = [S(1,:);S2d(2:end-1,:);S(end,:)]*[I;Rbar];
        b = Dmat\y;
        SS = sinfunc(n,L,xx);
        yp = (SS*[I;Rbar])*b;
        errvec(i,k) = errcompute(yp,yy);
        
        for j=1:N
            K(:,j) = cmatern(x,x(j),L,ep,beta);
            K2d(:,j) = cmatern(x,x(j),L,ep,beta,2);
            K_eval(:,j) = cmatern(xx,x(j),L,ep,beta);
        end
        K_solve = [K(1,:);K2d(2:end-1,:);K(end,:)];
        b = K_solve\y;
        yp = K_eval*b;
        errvecd(i,k) = errcompute(yp,yy);
        k = k+1;
    end
    i = i+1;
end

loglog(epvec,errvecd(1,:),'-bx')
hold on
loglog(epvec,errvecd(2,:),'-g+')
loglog(epvec,errvecd(3,:),'-r^')
loglog(epvec,errvec(1,:),'b','LineWidth',3)
loglog(epvec,errvec(2,:),'g','LineWidth',3)
loglog(epvec,errvec(3,:),'r','LineWidth',3)
hold off
xlabel('\epsilon')
ylabel('average error')
ptsstr=strcat(', x\in[',num2str(aa),',',num2str(bb),'],');
title(strcat(fstr,ptsstr,spacestr))
legend('N=10 (Direct)','N=20 (Direct)','N=40 (Direct)','N=10 (QR)','N=20 (QR)','N=40 (QR)', 'Location', 'SouthEast');

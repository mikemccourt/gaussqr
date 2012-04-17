% Tests RBF-QR to show that it can avoid the ill-conditioning
% associated with a small shape parameter
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

% May consider these functions which automatically satisfy the boundary
% conditions f(0)=f(L)=0
% yf = @(x) x.*(L-x);
% yf = @(x) x.*(L-x).*sqrt(x);

% Define the eigenfunctions and eigenvalues
sinfunc = @(n,L,x) sqrt(2/L)*sin(pi*x*n/L);
lamfunc = @(n,L,sigma,beta) ((pi*n/L).^2+sigma^2).^(-beta);

% Selecting some points in the domain, not including the boundary
aa = .1*L; bb = .9*L;
xx = pickpoints(aa,bb,NN);
yy = yf(xx);

% Set up the error vectors to store the results
errvec = zeros(length(Nvec),length(sigmavec));
errvecd = zeros(length(Nvec),length(sigmavecd));

% Choose the order of the basis functions
beta = 1;

i = 1;
for N=Nvec
    K_solve = zeros(N);
    K_eval = zeros(NN,N);
    [x,spacestr] = pickpoints(aa,bb,N,spaceopt);
    y = yf(x);
    I = eye(N);
    
    k = 1;
    for sigma=sigmavec
        M = ceil(8.5*N);
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
        errvec(i,k) = errcompute(yp,yy);
        
        for j=1:N
            K_solve(:,j) = sobfunc(x,x(j),L,sigma,beta);
            K_eval(:,j) = sobfunc(xx,x(j),L,sigma,beta);
        end
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

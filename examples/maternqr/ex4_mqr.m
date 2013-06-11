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

epvec =  logspace(-2,2,20);
Nvec = [10,20,40];
NN = 100;

spaceopt = 'even';

L = 1;
fstr = 'y(x)=exp(x)+x(1-e)-1';
yf = @(x) exp(x)+x*(1-exp(1))-1; % True solution
ff = @(x) exp(x); % Source

% Define the eigenvalues
lamfunc = mqr_solveprep();

% Selecting some points in the domain for computing the error
xx = pickpoints(0,L,NN);
yy = yf(xx);

% Set up the error vectors to store the results
errvec = zeros(length(Nvec),length(epvec));
errvecd = zeros(length(Nvec),length(epvec));

% Choose the order of the basis functions
beta = 2;

i = 1;
for N=Nvec
    K_solve = zeros(N);
    K_eval = zeros(NN,N);
    [x,spacestr] = pickpoints(0,L,N+2,spaceopt);
    x = x(2:end-1); % Dump 0 and L
    y = ff(x);
    I = eye(N);
    
    k = 1;
    for ep=epvec
        % Solve it directly first
        for j=1:N
            K_solve(:,j) = cmatern(x,x(j),L,ep,beta,2);
            K_eval(:,j) = cmatern(xx,x(j),L,ep,beta);
        end
        b = K_solve\y;
        yp = K_eval*b;
        errvecd(i,k) = errcompute(yp,yy);
        
        % Solve it with the QR method
        MQR = mqr_solveprep(x,L,ep,beta);
        Rbar = MQR.Rbar;
        Marr = MQR.Marr;
        
        S2d = mqr_phi(Marr,x,L,2);
        M_solve = S2d*[I;Rbar];
        M_eval = mqr_phi(Marr,xx,L);
        b = M_solve\y;
        yp = M_eval*([I;Rbar]*b);
        errvec(i,k) = errcompute(yp,yy);
        
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
ptsstr=strcat(', x\in[0,',num2str(L),'],');
title(strcat(fstr,ptsstr,spacestr))
legend('N=10 (Direct)','N=20 (Direct)','N=40 (Direct)','N=10 (QR)','N=20 (QR)','N=40 (QR)', 'Location', 'SouthEast');

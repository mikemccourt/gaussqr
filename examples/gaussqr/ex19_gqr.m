% ex18_gqr.m
% This should compute the likelihood function for a set of given data
% We are interested in looking at the relationship between this likelihood
% value and the error
global GAUSSQR_PARAMETERS

epvec = logspace(-2,1,31);

N = 15;
NN = 200;
x = pickpoints(-1,1,N,'cheb');
yf = @(x) x+1./(1+x.^2);
fstring = 'y(x) = x + 1/(1+x^2)';
yf = @(x) x.^3-3*x.^2+2*x+1;
fstring = 'y(x) = x^3-3x^2+2x+1';

y = yf(x);
xx = pickpoints(-1,1,NN);
yy = yf(xx);
alpha = 1;
lamratio = 1e-12;
pinvtol = 1e-11;

errvec = [];
mvec1  = [];
mvec2  = [];
cvec   = [];
derrvec = [];
dmvec   = [];
yPhi    = [];
yPsi    = [];
b       = [];
bPhi    = [];
detPhi1 = [];
detPsi  = [];
diffvec = [];


rbf = @(e,r) exp(-(e*r).^2);
DM = DistanceMatrix(x,x);
EM = DistanceMatrix(xx,x);

% Note that yPhi and yPsi are computed with a truncated SVD of Phi1 and Psi
% respectively.  The tolerance for this is pinvtol and can be set above.
k = 1;
for ep=epvec
    GQR = gqr_solve(x,y,ep,alpha,2*N+20);
    yp = gqr_eval(GQR,xx);
    errvec(k) = errcompute(yp,yy);
    
    Phi = gqr_phi(GQR,x);
    Phi1 = Phi(:,1:N);
    Phi2 = Phi(:,N+1:end);
    Psi = Phi1 + Phi2*GQR.Rbar;
    yPhi = Phi1\y;
    yPsi = Psi\y;
    beta = (1+(2*ep/alpha)^2)^.25;
    delta2 = alpha^2/2*(beta^2-1);
    ead = ep^2 + alpha^2 + delta2;
    lamvec = sqrt(alpha^2/ead)*(ep^2/ead).^(0:N-1)'; 
    laminv = 1./lamvec;
    lamsave = laminv.*(laminv/laminv(end)>lamratio);    
    Lambda1 = lamvec(1:N);
    Lambda2 = lamvec(N+1:end);
    b = Psi\y;
    bPhi = Phi1\Psi*b;
    
    %Mahalanobis Distance Calculation - Method One
    mahaldist1 = yPhi'*(lamsave.*yPsi);
    
    %Mahalanobis Distance Calculation - Method Two
    mahaldist2 = b'*(lamsave.*bPhi);
    
    %Mahalanobis Distance Vectors
    mvec1(k) = log(abs(mahaldist1));
    mvec2(k) = log(abs(mahaldist2));
    
    K = rbf(ep, DM);
    kbasis = rbf(ep,EM);
    warning off
    yp = kbasis*(K\y);    
    derrvec(k) = errcompute(yp,yy);
    
    %Condition vector of matrix K
    cvec(k) = cond(K);
    dmvec(k) = log(abs(y'*(K\y)));
    
    %Determinant of Phi1 and Psi
    detPhi1(k) = det(Phi1);
    detPsi(k) = det(Psi);
    
    %Difference between yPhi and yPsi
    diffvec(k) = sqrt(norm(abs(yPhi-yPsi)));
    
    warning on
    k = k + 1;
end

%Graph 1
loglog(epvec, exp(dmvec), 'g', 'linewidth', 3), hold on
loglog(epvec, cvec, 'b', 'linewidth', 3)
legend('direct norm','condition vector')
xlabel('\epsilon')
ylabel('Comparison of cond(K) and direct norm')
title(fstring), hold off
figure

%Graph 2
loglog(epvec, exp(dmvec), 'g', 'linewidth', 3), hold on
loglog(epvec, cvec, 'b', 'linewidth', 3)
loglog(epvec, exp(mvec1), 'm', 'linewidth', 3)
loglog(epvec, exp(mvec2), '--y', 'linewidth', 3)
legend('direct norm','condition vector', 'mvec1', 'mvec2')
xlabel('\epsilon')
ylabel('Comparison of cond(K), direct norm, mvec1, and mvec2')
title(fstring), hold off
figure

%Graph 3
semilogx(epvec, mvec1, 'm', 'linewidth', 3), hold on
semilogx(epvec, mvec2, '--y', 'linewidth', 3)
semilogx(epvec, dmvec, '--b', 'linewidth', 3)
legend('mvec1', 'mvec2', '-- direct')
xlabel('\epsilon')
ylabel('Comparison of Norms')
title(fstring), hold off
figure

%Graph 4
loglog(epvec, detPhi1, '--g', 'linewidth', 3), hold on
loglog(epvec, detPsi, '--b', 'linewidth', 3)
legend('--detPhi1', '--detPsi')
xlabel('\epsilon')
ylabel('Comparison of Determinants for \Phi_1 and \Psi')
title(fstring), hold off

%Graph 5
semilogx(epvec, diffvec, 'r', 'linewidth', 3), hold on
legend('diffvec')
xlabel('\epsilon')
ylabel('Difference between \y_Phi and \y_Psi')
title(fstring), hold off
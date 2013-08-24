% ex19_gqr.m
% This example will help us compare different methods of computation for
% the HS norm
global GAUSSQR_PARAMETERS

epvec = logspace(-2,1,31);

N = 800;
NN = 200;
x = pickpoints(-1,1,N,'cheb');
yf = @(x) x+1./(1+x.^2);
fstring = 'y(x) = x + 1/(1+x^2)';
% yf = @(x) x.^3-3*x.^2+2*x+1;
% fstring = 'y(x) = x^3-3x^2+2x+1';
fstring = sprintf('%s, N = %d',fstring,N);

y = yf(x);
xx = pickpoints(-1,1,NN);
yy = yf(xx);
alpha = 1;
lamratio = 1e-12;
pinvtol = 1e-11;

errvec = [];
mvec1  = [];
mvec2  = [];
mvec3  = [];
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
A       = [];
B       = [];

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
    yPsi = Psi\y;
    yPhi = Phi1\y;
    beta = (1+(2*ep/alpha)^2)^.25;
    delta2 = alpha^2/2*(beta^2-1);
    ead = ep^2 + alpha^2 + delta2;
    lamvec = sqrt(alpha^2/ead)*(ep^2/ead).^(0:N-1)';
    lamvec2 = sqrt(alpha^2/ead)*(ep^2/ead).^(N:size(GQR.Marr,2)-1)';
    laminv = 1./lamvec;
    lamsave = laminv.*(laminv/laminv(end)>lamratio);    
    Lambda1 = lamvec(1:N);
    Lambda2 = lamvec2;
    b = Psi\y;
    bPhi = Phi1\Psi*b;
%     B = (Phi2')/(Phi1')*diag(lamsave);
%     A = diag(lamsave) + B'*(diag(Lambda2))*B;
    
    
    %Mahalanobis Distance Calculation - Method One
    mahaldist1 = yPhi'*(lamsave.*yPsi);
    
    %Mahalanobis Distance Calculation - Method Two
    mahaldist2 = b'*(lamsave.*bPhi);
    
    %Mahalanobis Distance Calculation - Method Three (method two without
    %the correction term)z
    bvector = ((Lambda2.^(.5))'*(Phi2')/(Phi1')*(lamsave.*b))'*((Lambda2.^(.5))'*(Phi2')/(Phi1')*(lamsave.*b));
    mahaldist3 = b'*(lamsave.*b)+ bvector;
    
    %Print the MahalDist. using method 3
%     fprintf('%d\t%g\t%g\n',k,mahaldist3,bvector)

    %Mahalanobis Distance Calculation - Method Four (no cancellation; using
    %matrix A which is symmetric and pos. def.
    mahaldist4 = b'*(diag(lamsave) + ((Phi2')/(Phi1')*diag(lamsave))'*(diag(Lambda2))*(Phi2')/(Phi1')*diag(lamsave))*b;
    
    %Mahalanobis Distance Vectors
    mvec1(k) = log(abs(mahaldist1));
    mvec2(k) = log(abs(mahaldist2));
    mvec3(k) = log(abs(mahaldist3));
    mvec4(k) = log(abs(mahaldist4));
    
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
    
    %Comparison of Matrices
    absTol = 0.01;   % You choose this value to be what you want!
    relTol = 0.05;   % This one too!
    absError = Phi1(:)-Psi(:);
    relError = absError./Phi1(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = all( (abs(absError) < absTol) & (abs(relError) < relTol) );
    
    warning on
    k = k + 1;
end

% %Graph 1 - Comparison of condition vector and direct norm
% loglog(epvec, exp(dmvec), 'g', 'linewidth', 3), hold on
% loglog(epvec, cvec, 'b', 'linewidth', 3)
% legend('direct norm','condition vector')
% xlabel('\epsilon')
% ylabel('Comparison of cond(K) and direct norm')
% title(fstring), hold off
% figure
% 
% %Graph 2 - Comparison of Norms with the condition vector
% loglog(epvec, 10.^(dmvec), '--b', 'linewidth', 3), hold on
% loglog(epvec, cvec, 'g', 'linewidth', 3)
% loglog(epvec, 10.^(mvec1), 'm', 'linewidth', 3)
% loglog(epvec, 10.^(mvec2), '--y', 'linewidth', 3)
% loglog(epvec, 10.^(mvec3), '-.c', 'linewidth', 3)
% loglog(epvec, 10.^(mvec4), ':r', 'linewidth', 3)
% legend('direct norm','condition vector', 'mvec1', 'mvec2', 'mvec3', 'mvec4')
% xlabel('\epsilon')
% ylabel('Comparison of cond(K), direct norm, mvec1 through mvec4')
% title(fstring), hold off
% figure

%Graph 3 - Comparison of Norms without the condition vector
loglog(epvec, exp(mvec1), 'm', 'linewidth', 3), hold on
loglog(epvec, exp(mvec2), '--y', 'linewidth', 3)
loglog(epvec, exp(mvec3), '-.c', 'linewidth', 3)
loglog(epvec, exp(mvec4), ':r', 'linewidth', 3)
loglog(epvec, exp(dmvec), '--b', 'linewidth', 3)
legend('mvec1', 'mvec2', 'mvec3', 'mvec4', 'dmvec')
xlabel('\epsilon')
ylabel('Comparison of Norms')
title(fstring), hold off
% figure

% %Graph 4 - Comparison of Determinants for Phi1 and Psi
% loglog(epvec, detPhi1, '--g', 'linewidth', 3), hold on
% loglog(epvec, detPsi, '--b', 'linewidth', 3)
% legend('--det\Phi_1', '--det\Psi')
% xlabel('\epsilon')
% ylabel('Comparison of Determinants for \Phi_1 and \Psi')
% title(fstring), hold off
% figure
% 
% %Graph 5 - Difference between yPhi and yPsi
% loglog(epvec, diffvec, 'r', 'linewidth', 3), hold on
% legend('diffvec')
% xlabel('\epsilon')
% ylabel('Difference between y_\Phi and y_\Psi')
% title(fstring), hold off

beep

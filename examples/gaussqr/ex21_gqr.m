% ex21_gqr.m
% This example will involve an alternative method for computing the M. Dist
% We will utilize our Psi and b vectors in order to find y and from their
% determine our x and K
global GAUSSQR_PARAMETERS

epvec = logspace(-2,1,31);

N = 15;
NN = 200;
x = pickpoints(-1,1,N,'cheb');
fstring = 'b = ';
fstring = sprintf('%s, %d', fstring, b);

%Here is our vector b
b      = ones(N,1); %now B is defaulted to some vector
                    % of length N of all ones

alpha = 1;
lamratio = 1e-12;
pinvtol = 1e-11;

mvec1  = [];
mvec2  = [];
mvec3  = [];
cvec   = [];
dmvec   = [];
yPhi    = [];
yPsi    = [];
bPhi    = [];
detPhi1 = [];
detPsi  = [];
diffvec = [];
A       = [];
B       = [];
cvecPhi1= [];
cvecPsi = [];

rbf = @(e,r) exp(-(e*r).^2);
DM = DistanceMatrix(x,x);

% Note that yPhi and yPsi are computed with a truncated SVD of Phi1 and Psi
% respectively.  The tolerance for this is pinvtol and can be set above.
k = 1;
for ep=epvec
    GQR = gqr_solveprep(0,x,ep,alpha,2*N+20);
    
    Phi = gqr_phi(GQR,x);
    Phi1 = Phi(:,1:N);
    Phi2 = Phi(:,N+1:end);
    Psi = Phi1 + Phi2*GQR.Rbar;
    beta = (1+(2*ep/alpha)^2)^.25;
    delta2 = alpha^2/2*(beta^2-1);
    ead = ep^2 + alpha^2 + delta2;
    lamvec = sqrt(alpha^2/ead)*(ep^2/ead).^(0:N-1)';
    lamvec2 = sqrt(alpha^2/ead)*(ep^2/ead).^(N:size(GQR.Marr,2)-1)';
    laminv = 1./lamvec;
    lamsave = laminv.*(laminv/laminv(end)>lamratio);    
    Lambda1 = diag(lamvec(1:N));
    Lambda2 = diag(lamvec2);
    
    %approximate y value given b
    y   = Psi*b;
    c = Phi1'\Lambda1\b;
    
    %Mahalanobis Distance Calculation - Method One
    yPhi = Phi1\y;
    yPsi = Psi\y; %notice that this computation is actually equivalent
                % to the one for mahaldist2
    mahaldist1 = yPhi'*(lamsave.*yPsi);
    
    %Mahalanobis Distance Calculation - Method Two
    bPhi = Phi1\Psi*b;
    mahaldist2 = b'*(lamsave.*bPhi);
    
    %Mahalanobis Distance Calculation - Method Three (method two without
    %the correction term)z
    bvector = ((Lambda2.^(.5))'*(Phi2')/(Phi1')*(lamsave.*b))'*((Lambda2.^(.5))'*(Phi2')/(Phi1')*(lamsave.*b));
    mahaldist3 = b'*(lamsave.*b)+ bvector;
    
    %Mahalanobis Distance Calculation - Method Four (no cancellation; using
    %matrix A which is symmetric and pos. def.
    mahaldist4 = b'*(diag(lamsave) + ((Phi2')/(Phi1')*diag(lamsave))'*Lambda2*(Phi2')/(Phi1')*diag(lamsave))*b;
    
    %Mahalanobis Distance Vectors
    mvec1(k) = log(abs(mahaldist1));
    mvec2(k) = log(abs(mahaldist2));
    mvec3(k) = log(abs(mahaldist3));
    mvec4(k) = log(abs(mahaldist4));
    
    K = rbf(ep, DM);
    warning off
    
    %Condition vector of matrix K
    cvec(k) = cond(K);
    dmvec(k) = log(abs(y'*(K\y)));
    
    %Condition vectors for Phi1 and Psi
    cvecPhi1(k) = cond(Phi1);
    cvecPsi(k) = cond(Psi);
    
    %Determinant of Phi1 and Psi
    detPhi1(k) = det(Phi1);
    detPsi(k) = det(Psi);
    
    warning on
    k = k + 1;
end


%Graph 1 - Comparison of Norms with the condition vector
loglog(epvec, exp(dmvec), '--b', 'linewidth', 3), hold on
loglog(epvec, cvec, 'g', 'linewidth', 3)
loglog(epvec, exp(mvec1), 'm', 'linewidth', 3)
loglog(epvec, exp(mvec2), '--y', 'linewidth', 3)
loglog(epvec, exp(mvec3), '-.c', 'linewidth', 3)
loglog(epvec, exp(mvec4), ':r', 'linewidth', 3)
legend('direct norm','condition vector', 'mvec1', 'mvec2', 'mvec3', 'mvec4')
xlabel('\epsilon')
ylabel('Comparison of cond(K), direct norm, mvec1 through mvec4')
title(fstring), hold off
figure

%Graph 2 - Comparison of Norms without the condition vector
loglog(epvec, exp(mvec1), 'm', 'linewidth', 3), hold on
loglog(epvec, exp(mvec2), '--y', 'linewidth', 3)
loglog(epvec, exp(mvec3), '-.c', 'linewidth', 3)
loglog(epvec, exp(mvec4), ':r', 'linewidth', 3)
loglog(epvec, exp(dmvec), '--b', 'linewidth', 3)
legend('mvec1', 'mvec2', 'mvec3', 'mvec4', 'dmvec')
xlabel('\epsilon')
ylabel('Comparison of Norms')
title(fstring), hold off
figure

%Graph 6 - Comparison of Condition Numbers for Phi1, Psi, and K
loglog(epvec, cvecPhi1, '--g', 'linewidth', 3), hold on
loglog(epvec, cvecPsi, '--b', 'linewidth', 3)
loglog(epvec, cvec, '--r', 'linewidth', 3)
legend('cond(\Phi_1)','cond(\Psi)','cond(K)')
xlabel('\epsilon')
ylabel('Comparison of Condition Numbers for \Phi_1, \Psi, and K')
title(fstring), hold off

beep

%ex19a_gqr.m
%This example will help us continue our comparison that we saw in
%ex19_gqr.m, except that we will assign a value for epsilon and explore the
%plots for different values of N
global GAUSSQR_PARAMETERS

%Nvec = linspace(10, 200, 19);
Nvec = [10:1:40];

epsilon = 0.1; %we pick this
fstring = 'y(x) = x + 1/(1+x^2)'; %Function1
% fstring = 'y(x) = x^3-3x^2+2x+1'; %Function2
% fstring  = 'y(x) = 4tan(2x+6)'; %Function3
fstring = sprintf('%s, epsilon = %d',fstring,epsilon);


alpha = 1;
lamratio = 1e-12;
pinvtol = 1e-11;

mvec1  = [];
mvec2  = [];
mvec3  = [];
mvec4  = [];
cvec   = [];
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
cvecPhi1= [];
cvecPsi = [];

% Note that yPhi and yPsi are computed with a truncated SVD of Phi1 and Psi
% respectively.  The tolerance for this is pinvtol and can be set above.
k = 1;
for N=Nvec
    x = pickpoints(-1,1,N, 'cheb'); %how x depends on epsilon and not N...?
%    yf = @(x) x+1./(1+x.^2);        %Function1
    yf = @(x) x.^3-3*x.^2+2*x+1 + 0.001*cos(10*x);    %Function2
%     yf = @(x) 4*tan(2*x+6);         %Function3
    
    rbf=@(e,r) exp(-(e*r).^2);
    
    y = yf(x);
    GQR = gqr_solve(x,y,epsilon,alpha, N*2+20);
    
    DM = DistanceMatrix(x,x);
    
    Phi = gqr_phi(GQR,x);
    Phi1 = Phi(:,1:N);
    Phi2 = Phi(:,N+1:end);
    Psi = Phi1 + Phi2*GQR.Rbar;
    yPsi = Psi\y;
    yPhi = Phi1\y;
    beta = (1+(2*epsilon/alpha)^2)^.25;
    delta2 = alpha^2/2*(beta^2-1);
    ead = epsilon^2 + alpha^2 + delta2;
    lamvec = sqrt(alpha^2/ead)*(epsilon^2/ead).^(0:N-1)';
    lamvec2 = sqrt(alpha^2/ead)*(epsilon^2/ead).^(N:size(GQR.Marr,2)-1)';
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
    
    K = rbf(epsilon, DM);
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
    
    %Difference between yPhi and yPsi
    diffvec(k) = sqrt(norm(abs(yPhi-yPsi)));
    
    warning on
    k = k+1;
end

%Graph 1 - Comparison of Norms without the condition vector
semilogy(Nvec, exp(mvec1), 'm', 'linewidth', 3), hold on
semilogy(Nvec, exp(mvec2), '--y', 'linewidth', 3)
semilogy(Nvec, exp(mvec3), '-.c', 'linewidth', 3)
semilogy(Nvec, exp(mvec4), ':r', 'linewidth', 3)
semilogy(Nvec, exp(dmvec), '--b', 'linewidth', 3)
legend('mvec1', 'mvec2', 'mvec3', 'mvec4', 'dmvec')
xlabel('N')
ylabel('Comparison of Norms')
title(fstring), hold off
figure
semilogy(Nvec, cvecPhi1, 'm', 'linewidth', 3), hold on
semilogy(Nvec, cvecPsi, 'b', 'linewidth', 3), hold off

beep


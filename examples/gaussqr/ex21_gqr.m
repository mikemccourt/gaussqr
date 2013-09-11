% ex21_gqr.m
% This example will involve an alternative method for computing the M. Dist
% We will utilize our Psi and b vectors in order to find y and from their
% determine our x and K
global GAUSSQR_PARAMETERS

epvec = logspace(-2,1,31);

N = 15;
NN = 200;
x = pickpoints(-1,1,N,'cheb');
yf = @(x) x+1./(1+x.^2);
fstring = 'y(x) = x + 1/(1+x^2)';
% yf = @(x) x.^3-3*x.^2+2*x+1;
% fstring = 'y(x) = x^3-3x^2+2x+1';
% yf = @(x) 4*tan(2*x+6);
% fstring  = 'y(x) = 4tan(2x+6)';
fstring = sprintf('%s, N = %d',fstring,N);

%Here is our vector b
b      = ones(N,1); %now B is defaulted to some vector
                    % of length N of all ones

y = yf(x); %notice that this is negated below
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
EM = DistanceMatrix(xx,x);

% Note that yPhi and yPsi are computed with a truncated SVD of Phi1 and Psi
% respectively.  The tolerance for this is pinvtol and can be set above.
k = 1;
for ep=epvec
    %GQR = gqr_solve(x,y,ep,alpha,2*N+20);
    GQR = gqr_solveprep(0,x,ep,alpha)
    %yp = gqr_eval(GQR,xx);
    %errvec(k) = errcompute(yp,yy);
    
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
    pause
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
%     derrvec(k) = errcompute(yp,yy);
    
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

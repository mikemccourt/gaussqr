% ex17f_gqr.m
% This should compute the likelihood function for a set of given data points.
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION = 1;
GAUSSQR_PARAMETERS.ERROR_STYLE = 4;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

epvec = [logspace(-1,-.25,27),logspace(-.23,.8,80),logspace(.82,1,4)];
alpha = 1;
yf = @(x) cos(3*pi*x);

N = 20;
x = pickpoints(-1,1,N,'cheb');
x = pickpoints(-1,1,N);
y = yf(x);

NN = 100;
xx = pickpoints(-1,1,NN);
yy = yf(xx);

rbf = @(e,r) exp(-(e*r).^2);
DM_INT = DistanceMatrix(x,x);
DM_EVAL = DistanceMatrix(xx,x);

dirvec = zeros(size(epvec));
gqrvec = zeros(size(epvec));
errvec = zeros(size(epvec));

k = 1;
for ep=epvec
%     if ep<1
%         alpha = 1.4;
%     else
%         alpha = .6;
%     end
    % Solve the problem directly
    K = rbf(ep,DM_INT);
    S = svd(K);
    warning off
    dirvec(k) = N*log(y'*(K\y))+sum(log(S+eps));
    warning on
    
    % Solve the problem with HS-SVD
    GQR = gqr_solve(x,y,ep,alpha);
    errvec(k) = errcompute(gqr_eval(GQR,xx),yy);
    
    % Compute the log determinant
    Phi1 = GQR.stored_phi1;
    Phi2 = GQR.stored_phi2;
    S = svd(Phi1);
    logdetPhi = sum(log(S));
    
    Psi = Phi1 + Phi2*GQR.Rbar;
    S = svd(Psi);
    logdetPsi = sum(log(S));
       
    Lambda1 = GQR.eig(GQR.Marr(:,1:N));
    Lambda2 = GQR.eig(GQR.Marr(:,N+1:end));
     
    logdetK = logdetPsi + logdetPhi + sum(log(Lambda1));
    
    % Mahalanobis Distance
    laminv = 1./Lambda1;
    b = GQR.coef;
    mahaldist = (b'.*laminv)*b;
%     bvector = ((Lambda2.^(.5))'*(Phi2')/(Phi1')*(lamsave.*b));
%     mahaldist = b'*(laminv.*b)+ bvector'*bvector;
    gqrvec(k) = N*log(mahaldist) + logdetK;

    k = k + 1;
end

figure
[AX,H1,H2] = plotyy(epvec,[dirvec;gqrvec],epvec,errvec,'semilogx','loglog');
set(H1,'linewidth',3)
c = get(AX(1),'Children');
set(c(1),'color',[1 0 0])
set(c(1),'linestyle','--')
set(c(2),'color',[0 0 1])
set(H2,'linewidth',2)
set(H2,'color','k')
set(AX,{'ycolor'},{[.5 0 .5];[0 0 0]}) % Set axis color
legend('MLE Direct','MLE HS-SVD','Error','location','north')
xlabel('\epsilon')
ylabel('MLE function')
hold off
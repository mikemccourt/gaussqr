% ex17f_gqr.m
% This should compute the likelihood function for a set of given data
% points, but with a fixed coefficient vector b (as in Psi*b=y).
% Note that the vector y depends on epsilon and therefore changes in each 
% iteration. One probably needs to work with an adaptive choice of alpha
% here.
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION = 1;

epvec = logspace(-2,1,31);

N = 15;
x = pickpoints(-1,1,N,'cheb');
b = ones(size(x));
alpha = 1;
%lamratio = 1e-12;
lamratio = 0;
pinvtol = 1e-12;

detvec = [];
mvec   = [];
lvec   = [];
derrvec = [];
ddetvec = [];
dmvec   = [];
dlvec   = [];

rbf = @(e,r) exp(-(e*r).^2);
DM = DistanceMatrix(x,x);

k = 1;
for ep=epvec
    GQR = gqr_solveprep(0,x,ep);
    Phi1 = GQR.stored_phi1;
    Phi2 = GQR.stored_phi2;
    S = svd(Phi1);
    logdetPhi = sum(log(S));
    
    Psi = Phi1 + Phi2*GQR.Rbar;
    S = svd(Psi);
    logdetPsi = sum(log(S));
    
    beta = (1+(2*ep/alpha)^2)^.25;
    delta2 = alpha^2/2*(beta^2-1);
    ead = ep^2 + alpha^2 + delta2;
    Lambda1 = sqrt(alpha^2/ead)*(ep^2/ead).^(0:N-1)';
    Lambda2 = sqrt(alpha^2/ead)*(ep^2/ead).^(N:size(GQR.Marr,2)-1)';
     
    logdetK = logdetPsi + logdetPhi + sum(log(Lambda1));
    
    laminv = 1./Lambda1;
    lamsave = laminv.*(laminv/laminv(end)>lamratio);
    
    % Mahaldist
    bvector = ((Lambda2.^(.5))'*(Phi2')/(Phi1')*(lamsave.*b));
    mahaldist = b'*(lamsave.*b)+ bvector'*bvector;
    
    mvec(k) = log(abs(mahaldist));
    detvec(k) = 1/N*logdetK;
    lvec(k) = log(abs(mahaldist)) + 1/N*logdetK;

    y = Psi*b;
    figure, plot(x,y)
    A = rbf(ep,DM);
    warning off
    dmvec(k) = log(abs(y'*(A\y)));
    ddetvec(k) = 1/N*log(det(A));
    dlvec(k) =  dmvec(k) + ddetvec(k);
    warning on

    k = k + 1;
end

figure
semilogx(epvec,lvec,'r','linewidth',3), hold on
semilogx(epvec,dlvec,'--r','linewidth',3)
semilogx(epvec,mvec,'b','linewidth',3)
semilogx(epvec,dmvec,'--b','linewidth',3)
semilogx(epvec,detvec,'k','linewidth',3)
semilogx(epvec,ddetvec,'--k','linewidth',3)
legend('- log-like HS-SVD','- log-like direct','log(H_K-norm) HS','log(H_K-norm) direct','logdet(K)/N HS','logdet(K)/N direct')
xlabel('\epsilon')
ylabel('log-like function')
hold off
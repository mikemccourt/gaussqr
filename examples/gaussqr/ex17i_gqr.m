% ex17i_gqr.m
% This should compute the likelihood function for a set of given data
% points, but with a fixed coefficient vector b (as in Psi*b=y).
% Note that the vector y depends on epsilon and therefore changes in each 
% iteration. One probably needs to work with an adaptive choice of alpha
% here.
% We use both a range of epsilon and a range of N values (generalizes
% ex17f_gqr.m).
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION = 1;

epvec = logspace(-2,0,31);

Nvec = [5:5:50];
alpha = 1;
%lamratio = 1e-12;
lamratio = 0;

detvec = [];
mvec   = [];
lvec   = [];
derrvec = [];
ddetvec = [];
dmvec   = [];
dlvec   = [];

rbf = @(e,r) exp(-(e*r).^2);

l = 1;
for N=Nvec
    fprintf('N=%d \n',N)
    k = 1;
    x = pickpoints(-1,1,N,'cheb');
    DM = DistanceMatrix(x,x);
    %b = ones(N,1);
    b = zeros(N,1); b([1 2 5]) = 1;
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
    
        mvec(k,l) = log(abs(mahaldist));
        detvec(k,l) = 1/N*logdetK;
        lvec(k,l) = log(abs(mahaldist)) + 1/N*logdetK;

        y = Psi*b;
        A = rbf(ep,DM);
        warning off
        dmvec(k,l) = log(abs(y'*(A\y)));
        S = svd(A);
        ddetvec(k,l) = 1/N*sum(log(S));
        dlvec(k,l) =  dmvec(k,l) + ddetvec(k,l);
        warning on

        k = k + 1;
    end
    l = l + 1;
end

[X,Y] = meshgrid(epvec,Nvec);
surf(X,Y,lvec'), title('MLE')
set(gca,'XScale','log')
xlabel('\epsilon')
ylabel('N')
zlabel('MLE')
[i,j]=find(lvec==min(min(lvec)));
fprintf('MLE: optimal epsilon=%f, optimal N=%d\n',epvec(i),Nvec(j))

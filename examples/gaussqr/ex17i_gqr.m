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

epvec = logspace(-2,1,31);

%Nvec = [5:5:50];
Nvec = [10:10:200];
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
    b = zeros(N,1); b([1 2 5]) = 1;
%    b = zeros(N,1); b([1 2 10]) = 1;
%    b = ones(N,1);
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
        warning off
        bvector = ((Lambda2.^(.5))'*(Phi2')/(Phi1')*(lamsave.*b));
        warning on
        mahaldist = b'*(lamsave.*b)+ bvector'*bvector;
    
        lvec(k,l) = log(mahaldist) + 1/N*logdetK;

        y = Psi*b;
        A = rbf(ep,DM);
        warning off
        dmvec = log(abs(y'*(A\y)));
        S = svd(A);
        ddetvec = 1/N*sum(log(S));
        dlvec(k,l) =  dmvec + ddetvec;
        warning on

        k = k + 1;
    end
    l = l + 1;
end

[X,Y] = meshgrid(epvec,Nvec);
surf(X,Y,log10(exp(lvec'))), hold on, title('MLE - HS-SVD vs Direct')
surf(X,Y,log10(exp(dlvec')))
set(gca,'XScale','log')
xlabel('\epsilon')
ylabel('N')
zlabel('MLE')
for j=1:length(Nvec)
    i = find(lvec(:,j)==min(lvec(:,j)));
    scatter3(X(1,i),Y(j,1),log10(exp(lvec(i,j))),30,'fill','MarkerFaceColor','g')
end
[i,j]=find(lvec==min(min(lvec)));
fprintf('MLE HS-SVD: optimal epsilon=%f, optimal N=%d\n',epvec(i),Nvec(j))
scatter3(X(1,i),Y(j,1),log10(exp(lvec(i,j))),30,'fill','MarkerFaceColor','r')

% for j=1:length(Nvec)
%     i = find(dlvec(:,j)==min(dlvec(:,j)));
%     scatter3(X(1,i),Y(j,1),log10(exp(dlvec(i,j))),30,'fill','MarkerFaceColor','g')
% end
[i,j]=find(dlvec==min(min(dlvec)));
fprintf('MLE Direct: optimal epsilon=%f, optimal N=%d\n',epvec(i),Nvec(j))
%scatter3(X(1,i),Y(j,1),log10(exp(dlvec(i,j))),30,'fill','MarkerFaceColor','r')
hold off

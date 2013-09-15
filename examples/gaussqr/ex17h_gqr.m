% ex17h_gqr.m
% This should compute the likelihood function for a set of given data. We
% look at both a range of epsilon and N values.
% We are interested in looking at the relationship between this likelihood
% value and the error (generalizes ex17e_gqr.m).
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION = 1;

%epvec = logspace(-2,1,31);
epvec = logspace(-6,1,31);
%epvec = fliplr(epvec);

Nvec = [5:5:50];
NN = 200;
yf = @(x) x + 1./(1+x.^2);
fstring = 'y(x) = x + 1/(1+x^2)';
%yf = @(x) x.^3-3*x.^2+2*x+1 + 0.001*cos(10*x);
%fstring = 'y(x) = x^3-3x^2+2x+1';
%yf = @(x) x;% + 0.001*cos(10*x);
%fstring = 'y(x) = x';% + cos(10x)/1000';
%yf = @(x) 0.75*exp(-((9*(x+1)/2-2).^2)/4)+0.75*exp(-((9*(x+1)/2+1).^2/49))+0.5*exp(-((9*(x+1)/2-7).^2)/4)-0.2*exp(-((9*(x+1)/2-4).^2));
%fstring = '"Franke"';
%yf = @(x) tanh(9*(x-1))+1;
%fstring = 'tanh(9(x-1))+1';

% Why does this artificial "noise" make things work?

xx = pickpoints(-1,1,NN);
yy = yf(xx);
alpha = 1;
%lamratio = 1e-12;
lamratio = 0;

errvec = [];
detvec = [];
mvec   = [];
lvec   = [];
bvec   = [];
derrvec = [];
ddetvec = [];
dmvec   = [];
dlvec   = [];
lambdas = [];
lambdainv = [];
lambdasave = [];

rbf = @(e,r) exp(-(e*r).^2);

% Note that yPhi and yPsi are computed with a truncated SVD of Phi1 and Psi
% respectively.  The tolerance for this is pinvtol and can be set above.
l = 1;
for N=Nvec
    fprintf('N=%d \n',N)
    k = 1;
    x = pickpoints(-1,1,N,'cheb');
    y = yf(x);
    DM = DistanceMatrix(x,x);
    EM = DistanceMatrix(xx,x);
    for ep=epvec
        GQR = gqr_solve(x,y,ep,alpha);
        yp = gqr_eval(GQR,xx);
        errvec(k,l) = errcompute(yp,yy);
    
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
    
        % Mahaldist third version
        b = Psi\y;
        bvector = ((Lambda2.^(.5))'*(Phi2')/(Phi1')*(lamsave.*b));
        bvec(k,l) = bvector'*bvector;
        mahaldist = b'*(lamsave.*b) + bvec(k,l);
        mdist3(k,l) = mahaldist;
    
        mvec(k,l) = log(abs(mahaldist));
        detvec(k,l) = 1/N*logdetK;
        lvec(k,l) = log(abs(mahaldist)) + 1/N*logdetK;
    
        A = rbf(ep,DM);
        kbasis = rbf(ep,EM);
        warning off
        yp = kbasis*(A\y);
        derrvec(k,l) = errcompute(yp,yy);
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
surf(X,Y,errvec'), title('True solution')
set(gca,'XScale','log')
set(gca,'ZScale','log')
xlabel('\epsilon')
ylabel('N')
zlabel('error')
[i,j]=find(errvec==min(min(errvec)));
fprintf('True solution: optimal epsilon=%f, optimal N=%d\n',epvec(i),Nvec(j))

figure
[X,Y] = meshgrid(epvec,Nvec);
surf(X,Y,lvec'), title('MLE')
set(gca,'XScale','log')
xlabel('\epsilon')
ylabel('N')
zlabel('MLE')
[i,j]=find(lvec==min(min(lvec)));
fprintf('MLE: optimal epsilon=%f, optimal N=%d\n',epvec(i),Nvec(j))

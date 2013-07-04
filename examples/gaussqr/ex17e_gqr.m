% ex17e_gqr.m
% This should compute the likelihood function for a set of given data
% We are interested in looking at the relationship between this likelihood
% value and the error

epvec = logspace(-1,1,20);

N = 10;
NN = 200;
x = pickpoints(-1,1,N,'cheb');
yf = @(x) x+1./(1+x.^2);
%yf = @(x) x;
y = yf(x);
xx = pickpoints(-1,1,NN,'cheb');
yy = yf(xx);
alpha = 1;
lamratio = 1e-10;

errvec = [];
detvec = [];
mvec   = [];
lvec   = [];
k = 1;
for ep=epvec
    GQR = gqr_solve(x,y,ep,alpha);
    yp = gqr_eval(GQR,xx);
    errvec(k) = errcompute(yp,yy);
    
    Phi1 = GQR.stored_phi1;
    [U,S,V] = svd(Phi1);
    yPhi = V*((1./diag(S)).*(U'*y));
    logdetPhi = sum(log(diag(S)));
    
    Psi = Phi1 + GQR.stored_phi2*GQR.Rbar;
    [U,S,V] = svd(Psi);
    yPsi = V*((1./diag(S)).*(U'*y));
    logdetPsi = sum(log(diag(S)));

    beta = (1+(2*ep/alpha)^2)^.25;
    delta2 = alpha^2/2*(beta^2-1);
    ead = ep^2 + alpha^2 + delta2;
    lamvec = sqrt(alpha^2/ead)*(ep^2/ead).^(0:N-1)';
    
    logdetK = logdetPsi + logdetPhi + sum(log(lamvec));
    
    laminv = 1./lamvec;
    lamsave = laminv.*(laminv/laminv(end)>lamratio);
    
    mahaldist = yPhi'*(lamsave.*yPsi);
    
    mvec(k) = log(mahaldist);
    detvec(k) = 1/N*logdetK;
    lvec(k) = log(mahaldist) + 1/N*logdetK;
    k = k + 1;
end

loglog(epvec,errvec,'linewidth',3)
xlabel('\epsilon')
ylabel('error')
title('y(x) = x + 1/(1+x^2)')
figure
semilogx(epvec,[mvec;detvec;lvec],'linewidth',3)
legend('Mahal dist','determinant','log-like')
xlabel('\epsilon')
ylabel('log-like function ish')
title('y(x) = x + 1/(1+x^2)')
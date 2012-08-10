% Test studying the kriging variance
rbfsetup

rbf = @(e,r) exp(-(e*r).^2);

N = 10;
x = pickpoints(-1,1,N);
y = exp(x);

NN = 50;
xx = pickpoints(-1,1,NN);
yy = exp(xx);

xt = pickpoints(-.95,.95,N);

alpha = 1;

epvec = logspace(-1,.2,30);
kvvec = [];
errvec = [];
k = 1;
for ep=epvec
    GQR = gqr_solveprep(0,x,ep,alpha);
    phi = gqr_phi(GQR,x);
    phi1 = phi(:,1:N);
    psi = phi*[eye(N);GQR.Rbar];

    beta = (1+(2*ep/alpha)^2)^.25;
    delta2 = alpha^2/2*(beta^2-1);
    ead = ep^2 + alpha^2 + delta2;
    lamvec = sqrt(alpha^2/ead)*(ep^2/ead).^(0:N-1)';
    laminv = 1./lamvec;
    
    kvvec(k) = 0;
    for xv=xt'
        kx = rbf(ep,DistanceMatrix(x,xv));
        hkx = phi1\kx;
        skx = psi\kx;
        q = hkx'*(skx.*laminv);
        if q>rbf(ep,0)
            j = N-1;
            while q>rbf(ep,0) & j>0
                q = hkx'*(skx.*[laminv(1:j);zeros(N-j,1)]);
                j = j - 1
            end
            if q>rbf(ep,0)
                warning('unacceptable 1-q=%e, ep=%e, x=%e',1-q,ep,xv)
                q = rbf(ep,0);
            end
        end
%         A = rbf(ep,DistanceMatrix(x,x));
%         q = kx'*(A\kx);
        kvvec(k) = kvvec(k) + sqrt(rbf(ep,0)-q);
    end
    
    GQR = gqr_solve(x,y,ep,alpha);
    yp = gqr_eval(GQR,xx);
    errvec(k) = errcompute(yp,yy);
    
    k = k + 1;
end

loglog(epvec,[kvvec;errvec])
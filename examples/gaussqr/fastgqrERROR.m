% Look at the error when using various N and M combinations
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.NORM_TYPE = 2;
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;

f = @(x) besselj(0,6*(x+1));

N = 300;
NN = 1000;
x = pickpoints(-1,1,N,'halton');
xx = pickpoints(-1,1,NN);
y = f(x);
yy = f(xx);

Mvec = 4:2:36;
ep = .001;
alpha = 1;
opts.UT = true;
Mtol = 1e-16;

errorth_true = zeros(size(Mvec));
errorth_fast = zeros(size(Mvec));
errorth_fnew = zeros(size(Mvec));

errorth_Strue = zeros(size(Mvec));
errorth_Sfast = zeros(size(Mvec));
errorth_Sfnew = zeros(size(Mvec));
m = 1;
for M=Mvec
    I = eye(M);
    GQR = gqr_rsolve(x,y,ep,alpha,M);
    
    Marr = gqr_formMarr(M);
    phi = gqr_phi(Marr,x,ep,alpha);
    [Q,R] = qr(phi,0);
    errorth_true(m) = errcompute(Q'*Q,I);
    GQR.coef = R\(Q'*y);
    [yp,GQR] = gqr_eval(GQR,xx);
    errorth_Strue(m) = errcompute(yp,yy);
    
    [invU,Svec,Q] = computeQReig(M,x,ep,alpha);
    errorth_fast(m) = errcompute(Q'*Q,I);
    t = Q'*y;
    Mcut = min(find(abs(t)<Mtol));
    if isempty(Mcut)
        Mcut = M;
    end
    t = [t(1:Mcut);zeros(M-Mcut,1)];
    GQR.coef = invU*((1./Svec).*t);
    yp = gqr_eval(GQR,xx);
    errorth_Sfast(m) = errcompute(yp,yy);
    
    [invU,Svec,Q] = computeQReig_adjusted(M,x,ep,alpha);
    errorth_fnew(m) = errcompute(Q'*Q,I);
    t = Q'*y;
    Mcut = min(find(abs(t)<Mtol));
    if isempty(Mcut)
        Mcut = M;
    end
    t = [t(1:Mcut);zeros(M-Mcut,1)];
    GQR.coef = invU*((1./Svec).*t);
    yp = gqr_eval(GQR,xx);
    errorth_Sfnew(m) = errcompute(yp,yy);
    
    m = m+1;
end
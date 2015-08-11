% ex11c
% This example studies the viability and accuracy of the fast QR
% decomposition and least squares solver for the low-rank eigenfunction
% approximate basis

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
    GQR = gqr_solveprep(1,x,ep,alpha,M);
    
    [Q,R] = qr(gqr_phi(GQR,x),0);
    errorth_true(m) = errcompute(Q'*Q,I);
    GQR.coef = R\(Q'*y);
    [yp,GQR] = gqr_eval(GQR,xx);
    errorth_Strue(m) = errcompute(yp,yy);
    
    [invU,Svec,Q] = computeQReig(M,x,ep,alpha);
    errorth_fast(m) = errcompute(Q'*Q,I);
    t = Q'*y;
    Mcut = find(abs(t)<Mtol,1);
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
    Mcut = find(abs(t)<Mtol,1);
    if isempty(Mcut)
        Mcut = M;
    end
    t = [t(1:Mcut);zeros(M-Mcut,1)];
    GQR.coef = invU*((1./Svec).*t);
    yp = gqr_eval(GQR,xx);
    errorth_Sfnew(m) = errcompute(yp,yy);
    
    m = m+1;
end

subplot(1,2,1)
semilogy(Mvec,[errorth_true;errorth_fast;errorth_fnew])
xlabel('Eigenfunction series length')
ylabel('Error in orthogonality')
legend('True QR','Fast QR','Fast QR Adjusted')
subplot(1,2,2)
semilogy(Mvec,[errorth_Strue;errorth_Sfast;errorth_Sfnew])
xlabel('Eigenfunction series length')
ylabel('Error in solution')
legend('True QR','Fast QR','Fast QR Adjusted')
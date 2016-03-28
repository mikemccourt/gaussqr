% ex19
% This example is a demonstration of parameterization for a 2D tensor
% product kernel with one shape parameter per dimension, thus only 2 free
% parameters for the kernel.  This shows the stability afforded by the
% HS-SVD for small shape parameters.
% The initial demonstration, with significant comments explaining the
% mechanations, of the tensor product HS-SVD, are available in ex18.m.  We
% will avoid most of those comments here.

yf = @(x) cos(sqrt(sum(bsxfun(@times,x,[1,4]).^2,2)));

% Set up some data points at which to sample
N = 7;
x = pick2Dpoints(-1,1,N,'halton');
y = yf(x);

% Set up some evaluation points
NN = 10;
xx = pick2Dpoints(-1,1,NN);
yy = yf(xx);

% Parameters for the individual dimensions
% alpha is of little to no significance in interpolation
alpha = .5;
epvec = logspace(-4,-.3,25);

% The closed form of the Chebyshev kernel
K1d = @(e,x,z) alpha + (1-e)* ...
         (e*(1-e^2) - 2*e*bsxfun(@plus,x.^2,z.^2') + (1+3*e^2)*x*z')./ ...
         ((1-e^2)^2 + 4*e*(e*bsxfun(@plus,x.^2,z.^2')-(1+e^2)*x*z'));

% The Chebyshev kernel eigenfunctions and eigenvalues
phifunc = @(n,x) sqrt(2)*cos(acos(x)*n);
lambdafunc = @(e,n) (n==0).*(1-alpha) + (n>0).*(alpha*(1-e)/e*e.^n);

% Choose the eigenfunction buffer
% The actual value is irrelevant, just enough to create a Psi matrix
M = N + 50;
n1_test = 1:floor(M/5);
n2_test = 1:floor(M/5);
Marr_test = [kron(ones(size(n2_test)),n1_test);kron(n2_test,ones(size(n1_test)))];

errmat = zeros(length(epvec),length(epvec));
dirmat = zeros(length(epvec),length(epvec));
likmat = zeros(length(epvec),length(epvec));
dlimat = zeros(length(epvec),length(epvec));
k1 = 1;
for ep1=epvec
    k2 = 1;
    for ep2=epvec
        % Create the eigenvalue necessary ordering
        lamvec1_test = lambdafunc(ep1,Marr_test(1,:));
        lamvec2_test = lambdafunc(ep2,Marr_test(2,:));
        lamvec_test = lamvec1_test.*lamvec2_test;
        [lamsort,sorted_indices] = sort(lamvec_test,'descend');
        Marr = Marr_test(:,sorted_indices(1:M));
        
        % Create the HS-SVD basis
        lamvec = lamsort(1:M);
        lamvec1 = lamvec(1:N);
        lamvec2 = lamvec(N+1:M);
        Lambda1 = diag(lamvec1);
        Phi_K1 = phifunc(Marr(1,:),x(:,1));
        Phi_K2 = phifunc(Marr(2,:),x(:,2));
        Phi = Phi_K1.*Phi_K2;
        Phi1 = Phi(:,1:N);
        Phi2 = Phi(:,N+1:end);
        CbarT = (Phi2'/Phi1').*bsxfun(@rdivide,lamvec2',lamvec1);
        Psi = Phi1 + Phi2*CbarT;
        
        % Evaluate the interpolant at the relevant locations
        Phieval_K1 = phifunc(Marr(1,:),xx(:,1));
        Phieval_K2 = phifunc(Marr(2,:),xx(:,2));
        Phieval = Phieval_K1.*Phieval_K2;
        Psieval = Phieval*[eye(N);CbarT];
        b = Psi\y;
        yeval = Psieval*b;
        
        % Evaluate in that standard basis
        K = K1d(ep1,x(:,1),x(:,1)).*K1d(ep2,x(:,2),x(:,2));
        Keval = K1d(ep1,xx(:,1),x(:,1)).*K1d(ep2,xx(:,2),x(:,2));
        c = K\y;
        ydir = Keval*c;
        
        errmat(k1,k2) = errcompute(yeval,yy);
        dirmat(k1,k2) = errcompute(ydir,yy);
        
        % Compute profile likelihood in HS-SVD basis
        S = svd(Psi);
        logdetPsi = sum(log(S));
        S = svd(Phi1);
        logdetPhi = sum(log(S));
        logdetLam = sum(log(lamvec1));
        logdetK = logdetPsi + logdetPhi + logdetLam;
        L2P2P1L1invb = (lamvec2'.^-.5).*(CbarT*b);
        boundvec = b'*(b./lamvec1');
        correctionvec = L2P2P1L1invb'*L2P2P1L1invb;
        mahaldist = boundvec + correctionvec;
        likmat(k1,k2) = length(b)*log(mahaldist) + logdetK;
        
        % Compute profile likelihood in standard basis
        [U,S,V] = svd(K);
%         iii = find(abs(sum(U'*V-eye(length(U))))<1e-4);
%         Um = U(:,iii);
%         cm = c(iii);
%         Sm = S(iii,iii);
%         Spinv = [[inv(Sm),zeros(length(iii),length(S)-length(iii))];zeros(length(S)-length(iii),length(S))];
%         mahaldist = cm'*Um'*Spinv*Um*cm;
        mahaldist = abs(c'*y);
        logdetK = sum(log(diag(S)));
        dlimat(k1,k2) = length(c)*log(mahaldist) + logdetK;
        fprintf('%6.5f, %6.5f, %e, %e, %e\n',ep1, ep2, cond(K), likmat(k1,k2), dlimat(k1,k2))
        k2 = k2 + 1;
%         pause
    end
    k1 = k1 + 1;
end

f1 = figure;
[E1, E2] = meshgrid(epvec,epvec);
h = surf(E1, E2, errmat);
set(get(h,'parent'),'xscale','log');
set(get(h,'parent'),'yscale','log');
title('error in HS-SVD basis')

f2 = figure;
h = surf(E1, E2, dirmat);
set(get(h,'parent'),'xscale','log');
set(get(h,'parent'),'yscale','log');
title('error in standard basis')

f3 = figure;
h = surf(E1, E2, likmat);
set(get(h,'parent'),'xscale','log');
set(get(h,'parent'),'yscale','log');
title('likelihood in HS-SVD basis')

f4 = figure;
h = surf(E1, E2, dlimat);
set(get(h,'parent'),'xscale','log');
set(get(h,'parent'),'yscale','log');
title('likelihood in standard basis')
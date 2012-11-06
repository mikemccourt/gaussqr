N = 10;
x = pickpoints(0,1,N,'cheb');
yf = @(x) 30*x.^2.*(1-x).^2.*sin(2*pi*x).^4;
y = yf(x);        
% Choice of kernel
ep = 10;
beta = 5;   % DO NOT use with beta < 3 !!
% Eigenvalues and eigenfunctions
lamfunc = @(n,ep,beta) ((pi*n).^2+ep^2).^(-beta);
phifunc = @(n,x) sqrt(2)*sin(pi*x*n);
% Truncation length
M = floor(1/pi*sqrt(eps^(-1/beta)*(pi^2+ep^2)-ep^2));
n = 1:M;
% Fill the eigenfunction matrix
Phi = phifunc(n,x);
% Compute the QR decomposition of the short, fat matrix Phi
[Q,R] = qr(Phi);
R1 = R(:,1:N);
R2 = R(:,N+1:end);
% Apply inv(R1)*R2
opts.UT = true;
Rhat = linsolve(R1,R2,opts);
% Compute the eigenvalue matrix
lambda = lamfunc(n,ep,beta);
D1 = repmat(lambda(1:N),M-N,1);
D2 = repmat(lambda(N+1:end)',1,N);
% Form the Rbar correction matrix
Rbar = D2.*Rhat'./D1;
% Solve for the interpolant coefficients
I = eye(N);
coef = linsolve(Phi*[I;Rbar],y);
% Evaluate the approximation
xx = pickpoints(0,1,500,'even');
phiEval1 = phifunc(1:N,xx);
phiEval2 = phifunc(N+1:M,xx);
yy = phiEval1*coef+phiEval2*Rbar*coef;
% Plot solution
plot(xx,yy)
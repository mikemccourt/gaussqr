% MaternQRDemo
% Computes and plot a stable and unstable compact Matern interpolant
ep = 1; beta = 7;
G14 = @(x) .25^(-28)*max(x-.25,0).^14.*max(.75-x,0).^14;
N = 21;
x = pickpoints(0,1,N); x = x(2:N-1); y = G14(x); N = N-2;
phifunc = @(n,x) sqrt(2)*sin(pi*x*n);
M = N + floor(1/pi*sqrt(eps^(-1/beta)*(pi^2+ep^2)-ep^2));
% Create a vector of the desired eigenvalues
n = 1:M;
Lambda = diag(((pi*n).^2+ep^2).^(-beta));
Phi = phifunc(n,x);
% First solve in the standard basis
K = Phi*Lambda*Phi';
c = K\y;
% Now solve with the RBF-QR technique
[Q,R] = qr(Phi);
R1 = R(:,1:N); R2 = R(:,N+1:end);
Rhat = R1\R2;
Lambda1 = Lambda(1:N,1:N); Lambda2 = Lambda(N+1:M,N+1:M);
Rbar = Lambda2*Rhat'/Lambda1;
Psi = Phi*[eye(N);Rbar];
b = Psi\y;
% Evaluate the results and plot them
x_eval = linspace(0,1,100)';
Phi_eval = phifunc(n,x_eval);
y_standard = Phi_eval*Lambda*Phi'*c;
y_RBFQR = Phi_eval*[eye(N);Rbar]*b;
plot(x_eval,y_standard,'linewidth',2),hold on
plot(x_eval,y_RBFQR,'g',x_eval,G14(x_eval),':r','linewidth',3)
legend('Standard basis','MaternQR','True solution'),hold off
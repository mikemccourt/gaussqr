% ChebQRDemo
% Computes and plots a stable and unstable Chebyshev kernel interpolant
alpha = .5; beta = .1;
testfun = @(x) 4*sin(3*x)./exp(3*x/2);
N = 51;
x = -cos(pi*(0:N-1)/(N-1))';
y = testfun(x);
phifunc = @(n,x) sqrt(2)*cos(acos(x)*n);
M = N + ceil(log(eps)/log(beta));
% Create a vector of the desired eigenvalues
n = 1:M;
Lambda = diag([1-alpha alpha*(1-beta)/beta*beta.^n]);
Phi = [ones(N,1) phifunc(n,x)];
% First solve in the standard basis
K = Phi*Lambda*Phi';
c = K\y;
% Now solve with the RBF-QR technique
[Q,R] = qr(Phi);
R1 = R(:,1:N); R2 = R(:,N+1:end);
Rhat = R1\R2;
Lambda1 = Lambda(1:N,1:N); Lambda2 = Lambda(N+1:end,N+1:end);
Rbar = Lambda2*Rhat'/Lambda1;
Psi = Phi*[eye(N);Rbar];
b = Psi\y;
% Evaluate the results and plot them
x_eval = linspace(-1,1,100)';
Phi_eval = [ones(100,1) phifunc(n,x_eval)];
y_standard = Phi_eval*Lambda*Phi'*c;
y_RBFQR = Phi_eval*[eye(N);Rbar]*b;
y_eval = testfun(x_eval);
plot(x_eval,y_standard,'linewidth',2),hold on
plot(x_eval,y_RBFQR,'g',x_eval,y_eval,':r','linewidth',3)
legend('Standard basis','ChebQR','True solution'),hold off
figure
err_standard = abs(y_standard-y_eval);
err_RBFQR = abs(y_RBFQR-y_eval)+eps;    % There seem to be a few places with 
                                        % zero error, so to make this look 
                                        % nicer a little shift by eps 
semilogy(x_eval,err_standard,'linewidth',2),hold on
semilogy(x_eval,err_RBFQR,'g','linewidth',3)
legend('Standard basis','ChebQR'),hold off

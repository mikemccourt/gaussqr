% MaternQRLimit
% Studies the stable and unstable ep->0 limit for compact Matern
% Calls on: errcompute, pickpoints
N_vals = [15,30,45]; ep_vals = logspace(0,2,40);
beta = 7;
phifunc = @(n,x) sqrt(2)*sin(pi*x*n);
G14 = @(x) .25^(-28)*max(x-.25,0).^14.*max(.75-x,0).^14;
x_eval = linspace(0,1,200)'; y_eval = G14(x_eval);
err_RBFQR = zeros(length(N_vals),length(ep_vals));
err_standard = zeros(length(N_vals),length(ep_vals));
j = 1;
for N=N_vals
    x = pickpoints(0,1,N+2); x = x(2:N-1); y = G14(x);
    k = 1;
    for ep=ep_vals
        M = N + floor(1/pi*sqrt(eps^(-1/beta)*(pi^2+ep^2)-ep^2));
        n = 1:M;
        Lambda = diag(((pi*n).^2+ep^2).^(-beta));
        Phi = phifunc(n,x);
        Phi_eval = phifunc(n,x_eval);
        % Solve the problem directly
        K = Phi*Lambda*Phi'; c = K\y;
        y_standard = Phi_eval*Lambda*Phi'*c;
        err_standard(j,k) = errcompute(y_standard,y_eval);
        % Solve the problem with RBF-QR
        [Q,R] = qr(Phi);
        R1 = R(:,1:N); R2 = R(:,N+1:end);
        Rhat = R1\R2;
        lvec = ((pi*n).^2+ep^2).^(-beta);
        Lambda21 = repmat(lvec(N+1:M)',1,N)./repmat(lvec(1:N),M-N,1);
        RbarT = Lambda21.*Rhat';
        Lambda1 = Lambda(1:N,1:N); Lambda2 = Lambda(N+1:M,N+1:M);
        Rbar = Lambda2*Rhat'/Lambda1;
        lambda = lvec;
D1 = repmat(lambda(1:N),M-N,1);
D2 = repmat(lambda(N+1:end)',1,N);

% Form the Rbar matrix
Rbar = D2.*Rhat'./D1;
        fprintf('%d\t%d\t%g\n',j,k,norm(Rbar-RbarT)/norm(RbarT))
        Psi = Phi*[eye(N);Rbar]; b = Psi\y;
        y_RBFQR = Phi_eval*[eye(N);Rbar]*b;
        y_RBFQR = mqr_eval(mqr_solve(x,y,1,ep,beta),x_eval);
        err_RBFQR(j,k) = errcompute(y_RBFQR,y_eval);
        k = k + 1;
    end
    j = j + 1;
end
loglog(ep_vals,err_standard,'--','linewidth',3),hold on
loglog(ep_vals,err_RBFQR,'linewidth',3),hold off
legend('N=15','N=30','N=45')
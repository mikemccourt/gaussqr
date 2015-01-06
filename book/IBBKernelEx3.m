% IBBkernelEx3.m
% To illustrate convergence rates as a function of beta
% Calls on: errcompute, pickpoints, IBBkernel (for beta=1,2)
ep = 1; betavec = 1:5; Nvec = 2.^(3:8); Neval = 397;
f = @(x) 5.0887^6.*(max(0,x-.0567)).^6.*(max(0,(1-.0567)-x)).^6;
x_eval = linspace(0,1,Neval)'; y_eval = f(x_eval);
errvec = zeros(length(betavec),length(Nvec));
errorForBetas = zeros(length(betavec),Neval);
k = 1; j = 1;
for N=Nvec
    x = pickpoints(0,1,N+2); x = x(2:end-1); y = f(x);
    phifunc = @(n,x) sqrt(2)*sin(pi*x*n);
    j = 1;
    for beta=betavec
        if beta==1 || beta==2 % Work with the kernel form
            K_solve = zeros(N); K_eval = zeros(Neval,N);
            for n=1:N 
                K_solve(:,n) = IBBkernel(x,x(n),1,ep,beta);
                K_eval(:,n) = IBBkernel(x_eval,x(n),1,ep,beta);
            end
            c = K_solve\y;
            yp = K_eval*c; 
        else % Work with the series form
            M = ceil(1/pi*sqrt(eps^(-1/beta)*(pi^2+ep^2)-ep^2)); 
            M = max(M,N+1); 
            n = 1:M; 
            Lambda = diag(((pi*n).^2+ep^2).^(-beta));
            Phi = phifunc(n,x);
            [Q,R] = qr(Phi);
            R1 = R(:,1:N); R2 = R(:,N+1:end);
            Rhat = R1\R2;
            Lambda1 = Lambda(1:N,1:N);
            Lambda2 = Lambda(N+1:M,N+1:M);
            Rbar = Lambda2*Rhat'/Lambda1;
            Psi = Phi*[eye(N);Rbar];
            b = Psi\y;
            Phi_eval = phifunc(n,x_eval);
            yp = Phi_eval*[eye(N);Rbar]*b; 
        end
        errvec(j,k) = errcompute(yp,y_eval);
        errorForBetas(beta,:) = y_eval-yp;
        j = j + 1;
    end
    k = k + 1;
end
% Plot RMS error as beta and N vary
figure, loglog(Nvec,errvec,'linewidth',2); 
% Plot individual error
figure, title('Error as \beta varies'); 
for beta = betavec;
    errorFig = subplot(ceil(length(betavec)/2),2,beta);
    semilogy(x_eval,abs(errorForBetas(beta,:)));
    ylim([1e-18 1e-3]); set(gca,'YTick',[1e-15 1e-10 1e-5])
end 

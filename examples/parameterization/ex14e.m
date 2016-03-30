% ex14e.m
% This computes various parametrization criteria for a range of b and 
% fixed input data using Chebyshev kernels and the HS-SVD for stability in 
% smaller b.
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 4;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

a = 0.4;
%epvec = logspace(-2,-.01,80);
epvec = [logspace(-2,-.5,20),logspace(-.5,-.01,60)];
%yf = @(x) cos(3*pi*x);
%yf = @(x) 4*sin(4*pi*x)./exp(3*x/2);
yf = @(x) x.*log(1+x.^2);

N = 24;
%x = pickpoints(-1,1,N);
x = pickpoints(-1,1,N,'cheb');
y = yf(x);

NN = 50;
xx = pickpoints(-1,1,NN);
yy = yf(xx);

% Define eigenfunctions and eigenvalues of the rational Chebyshev kernel
% with geometrically decaying eigenvalues
phifunc = @(n,x) sqrt(2)*cos(acos(x)*n);

dirMPLEvec = zeros(size(epvec));
dirKVvec = zeros(size(epvec));
dirDETvec = zeros(size(epvec));
gqrMPLEvec = zeros(size(epvec));
gqrKVvec = zeros(size(epvec));
gqrKVdetvec = zeros(size(epvec));
gqrDETvec = zeros(size(epvec));
gqrDETdetvec = zeros(size(epvec));
errvec = zeros(size(epvec));
boundvec = zeros(size(epvec));
correctionvec = zeros(size(epvec));
denvec = zeros(size(epvec));
nummat = zeros(length(NN),length(epvec));

% Decide if we want to use the max-norm or 1-norm
KV_max_norm = 1;

k = 1;
h_waitbar = waitbar(0,'Initializing');
for ep=epvec
    % Solve the problem directly
    M = N + ceil(log(eps)/log(ep));
    n = 1:M;
    Lambda = diag([1-a a*(1-ep)*ep.^(n-1)]);
    Phi = [ones(N,1) phifunc(n,x)];
    % First solve in the standard basis
    K = Phi*Lambda*Phi';
    [U,S,V] = svd(K);
    Mdist = y'*(V*diag(0*(diag(S)<1e-14)+(1./diag(S)).*(diag(S)>1e-14))*U'*y);
    logdetK = sum(log(diag(S)+eps));
    Phi_eval = [ones(NN,1) phifunc(n,xx)];
    Keval = Phi_eval*Lambda*Phi';
    Kx0 = diag(Phi_eval*Lambda*Phi_eval');
    warning off
    pvals = Kx0 - sum((Keval/K).*Keval,2);
    warning on
    if KV_max_norm
        PF = max(pvals);
    else
        PF = sum(abs(pvals));
    end
    dirMPLEvec(k) = N*log(Mdist) + logdetK;
    dirKVvec(k) = log(Mdist) + log(PF);
    dirDETvec(k) = N*(dirMPLEvec(k) + dirKVvec(k))/(N+1); % multiply by N/(N+1) to have same scale

    % Now solve with the HS-SVD technique
    [Q,R] = qr(Phi);
    R1 = R(:,1:N); R2 = R(:,N+1:end);
    Rhat = R1\R2;
    Lambda1 = Lambda(1:N,1:N); Lambda2 = Lambda(N+1:end,N+1:end);
    lamvec1 = diag(Lambda1); lamvec2 = diag(Lambda2);
    Rbar = Rhat'.*bsxfun(@rdivide,lamvec2,lamvec1');
    Psi = Phi*[eye(N);Rbar];
    b = Psi\y;
    Psi_eval = Phi_eval*[eye(N);Rbar]; 
    y_eval = Psi_eval*b;
    errvec(k) = errcompute(y_eval,yy);
    
    Lambda1vec = diag(Lambda1); Lambda2vec = diag(Lambda2);
    
    % Compute the log determinant
    Phi1 = Phi(:,1:N);
    Phi2 = Phi(:,N+1:end);
    S = svd(Phi1);
    logdetPhi = sum(log(S));

    Psi = Phi1 + Phi2*Rbar;
    S = svd(Psi);
    logdetPsi = sum(log(S));
     
    logdetK = logdetPsi + logdetPhi + sum(log(Lambda1vec));
    
    % Mahalanobis Distance
    L2P2P1L1invb = (Lambda2vec.^-.5).*(Rbar*b);
    boundvec(k) = b'*(b./Lambda1vec);
    correctionvec(k) = L2P2P1L1invb'*L2P2P1L1invb;
    mahaldist = boundvec(k) + correctionvec(k);

    % power function
    pvals = Kx0 - sum((Psi_eval/Psi).*Keval,2);
    if KV_max_norm
        PF = max(pvals);
    else
        PF = sum(abs(pvals));
    end

    gqrMPLEvec(k) = N*log(mahaldist) + logdetK;
    gqrKVvec(k) = log(mahaldist) + log(PF);

    % Using determinant identity to compute power function without cancelation
    % First the denominator
    denvec(k) = logdetK;
    % Then the numerator
    for m=1:NN
        xp = [x;xx(m)];
        Np = N + 1;
        Mp = Np + ceil(log(eps)/log(ep));
        np = 1:Mp;
        Lambda = diag([1-a a*(1-ep)*ep.^(np-1)]);
        Phi = [ones(Np,1) phifunc(np,xp)];
        Phi1 = Phi(:,1:Np);
        Phi2 = Phi(:,Np+1:end);
        [Q,R] = qr(Phi);
        R1 = R(:,1:Np); R2 = R(:,Np+1:end);
        Rhat = R1\R2;
        Lambda1 = Lambda(1:Np,1:Np); Lambda2 = Lambda(Np+1:end,Np+1:end);
        lamvec1 = diag(Lambda1); lamvec2 = diag(Lambda2);
        Rbar = Rhat'.*bsxfun(@rdivide,lamvec2,lamvec1');
        S = svd(Phi1);
        logdetPhi = sum(log(S));
        Psi = Phi1 + Phi2*Rbar;
        S = svd(Psi);
        logdetPsi = sum(log(S));
        Lambda1vec = diag(Lambda1);
        nummat(m,k) = logdetPsi + logdetPhi + sum(log(Lambda1vec));
    end
    logpvals = nummat(:,k) - denvec(k);
    if KV_max_norm
        norm_log_Pf = max(logpvals);
        norm_log_det_Ktilde = max(nummat(:,k));
    else
        norm_log_Pf = log(sum(exp(logpvals))); % 1-norm of vector of values of power function squared
        norm_log_det_Ktilde = log(sum(exp(nummat(:,k))));
    end

    gqrKVdetvec(k) = log(mahaldist) + norm_log_Pf; % multiply by N to have same scale
    gqrDETvec(k) = N*(gqrMPLEvec(k) + gqrKVdetvec(k))/(N+1); % multiply by N/(N+1) to have same scale
    gqrDETdetvec(k) = N*((N+1)*log(mahaldist) + norm_log_det_Ktilde)/(N+1); % multiply by N/(N+1) to have same scale
    
    progress = floor(100*k/length(epvec))/100;
    waitbar(progress,h_waitbar,sprintf('Computing, ep=%g',ep))
    k = k + 1;
end

waitbar(100,h_waitbar,sprintf('Plotting'))
h_mle = figure;
[AX,H1,H2] = plotyy(epvec,[dirMPLEvec;gqrMPLEvec;N*dirKVvec;N*gqrKVvec;N*gqrKVdetvec;dirDETvec;gqrDETvec;gqrDETdetvec],epvec,errvec,'semilogx','loglog');
set(H1,'linewidth',3)
c = get(AX(1),'Children');
set(c(1),'color',[0.9290 0.6940 0.1250])
set(c(1),'linestyle',':')
set(c(2),'color',[0.9290 0.6940 0.1250])
set(c(2),'linestyle','-')
set(c(3),'color',[0.9290 0.6940 0.1250])
set(c(3),'linestyle','--')
set(c(4),'color',[0.8500 0.3250 0.0980])
set(c(4),'linestyle','-')
set(c(5),'color',[0.8500 0.3250 0.0980])
set(c(5),'linestyle',':')
set(c(6),'color',[0.8500 0.3250 0.0980])
set(c(6),'linestyle','--')
set(c(7),'color',[0 0.4470 0.7410])
set(c(7),'linestyle','-')
set(c(8),'color',[0 0.4470 0.7410])
set(c(8),'linestyle','--')
set(H2,'linewidth',2)
set(H2,'color','k')
set(AX,{'ycolor'},{[.5 0 .5];[0 0 0]}) % Set axis color
set(AX(1),'ylim',[-1500,1500])
set(AX(2),'ylim',[1e-18,1e-2])
set(AX(2),'ytick',[1e-9,1e-6,1e-3])
legend([H1;H2],'MLE direct','MLE HS-SVD','KV direct (x N)','KV HS-SVD (x N)','KV HS-SVDdet (x N)','DET direct','DET HS-SVD','DET HS-SVDdet','Error','location','north')
xlabel('\epsilon')
set(get(AX(1),'xlabel'),'Position',[0.3,-1560,0])
ylabel('Likelihood function')
set(get(AX(1),'ylabel'),'Position',[.08,0,17])
set(get(AX(2),'ylabel'),'String','Relative error')
set(get(AX(2),'ylabel'),'Position',[1.1 3e-014 0])
hold off

close(h_waitbar)
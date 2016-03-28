% ex14b.m
% This computes the likelihood for a range of epsilon and fixed input data
% using the HS-SVD for stability in smaller epsilon.
% It also analyzes the effectiveness of the approximate bound, rather than
% the full computation.
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION = 1;
GAUSSQR_PARAMETERS.ERROR_STYLE = 4;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

epvec = [logspace(-1,-.25,27),logspace(-.23,.8,80),logspace(.82,1,4)];
alpha = 2;
yf = @(x) cos(3*pi*x);

N = 24;
x = pickpoints(-1,1,N);
% x = pickpoints(-1,1,N,'cheb');
y = yf(x);

NN = 50;
xx = pickpoints(-1,1,NN);
yy = yf(xx);

rbf = @(e,r) exp(-(e*r).^2);
DM_INT = DistanceMatrix(x,x);
DM_EVAL = DistanceMatrix(xx,x);

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
    K = rbf(ep,DM_INT);
    [U,S,V] = svd(K);
    Mdist = y'*(V*diag(0*(diag(S)<1e-14)+(1./diag(S)).*(diag(S)>1e-14))*U'*y);
    logdetK = sum(log(diag(S)+eps));
    Keval = rbf(ep,DM_EVAL);
    warning off
    pvals = 1 - sum((Keval/K).*Keval,2);
    warning on
    if KV_max_norm
        PF = max(pvals);
    else
        PF = sum(abs(pvals));
    end
    dirMPLEvec(k) = N*log(Mdist) + logdetK;
    dirKVvec(k) = log(Mdist) + log(PF); % multiply by N to have same scale
    dirDETvec(k) = N*(dirMPLEvec(k) + dirKVvec(k))/(N+1); % multiply by N/(N+1) to have same scale
    
    % Solve the problem with HS-SVD
    GQR = gqr_solve(x,y,ep,alpha);
    errvec(k) = errcompute(gqr_eval(GQR,xx),yy);
    
    % Compute the log determinant
    Phi1 = GQR.stored_phi1;
    Phi2 = GQR.stored_phi2;
    Rbar = GQR.Rbar;
    S = svd(Phi1);
    logdetPhi = sum(log(S));
    
    Psi = Phi1 + Phi2*Rbar;
    S = svd(Psi);
    logdetPsi = sum(log(S));
       
    Lambda1 = GQR.eig(GQR.Marr(:,1:N))';
    Lambda2 = GQR.eig(GQR.Marr(:,N+1:end))';
     
    logdetK = logdetPsi + logdetPhi + sum(log(Lambda1));
    
    % Mahalanobis Distance
    b = GQR.coef;
    L2P2P1L1invb = (Lambda2.^-.5).*(Rbar*b);
    boundvec(k) = b'*(b./Lambda1);
    correctionvec(k) = L2P2P1L1invb'*L2P2P1L1invb;
    mahaldist = boundvec(k) + correctionvec(k);

    % power function
    Phieval = gqr_phi(GQR,xx);
    Psieval = Phieval*[eye(N);Rbar];
    pvals = 1 - sum((Psieval/Psi).*Keval,2);
    if KV_max_norm
        PF = max(pvals);
    else
        PF = sum(abs(pvals)); % 1-norm of vector of values of power function squared
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
        GQR = gqr_solveprep(0,xp,ep,alpha);
        Phi1 = GQR.stored_phi1;
        Phi2 = GQR.stored_phi2;
        Rbar = GQR.Rbar;
        S = svd(Phi1);
        logdetPhi = sum(log(S));
        Psi = Phi1 + Phi2*Rbar;
        S = svd(Psi);
        logdetPsi = sum(log(S));
        Lambda1 = GQR.eig(GQR.Marr(:,1:Np))';
        nummat(m,k) = logdetPsi + logdetPhi + sum(log(Lambda1));
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
set(AX(1),'ylim',[-1000,1500])
set(AX(2),'ylim',[1e-18,1e-2])
set(AX(2),'ytick',[1e-9,1e-6,1e-3])
legend([H1;H2],'MLE direct','MLE HS-SVD','KV direct (x N)','KV HS-SVD (x N)','KV HS-SVDdet (x N)','DET direct','DET HS-SVD','DET HS-SVDdet','Error','location','north')
xlabel('\epsilon')
set(get(AX(1),'xlabel'),'Position',[2,-1060,0])
ylabel('Likelihood function')
set(get(AX(1),'ylabel'),'Position',[.07,500,17])
set(get(AX(2),'ylabel'),'String','Relative error')
set(get(AX(2),'ylabel'),'Position',[12 3e-014 0])
hold off

close(h_waitbar)
% ex19b
% This is a 2D example for stable parametrization using both the MLE and
% the joint metric.
% We use the Gaussian in 2D and try to optimize for a single ep value.

yf = @(x) cos(sqrt(sum(bsxfun(@times,x,[.5,1]).^2,2)));

% Set up some data points at which to sample
N = 9;
x = pick2Dpoints(-1,1,N,'halton');
y = yf(x);

% Set up some evaluation points
NN = 10;
xx = pick2Dpoints(-1,1,NN);
yy = yf(xx);
X = reshape(xx(:,1),[NN,NN]);
Y = reshape(xx(:,2),[NN,NN]);

% Parameters for the individual dimensions
alpha = 1;
epvec = [logspace(-2,-1,8),logspace(-1,0,50),logspace(0,1,8)];

% The closed form of the radial Gaussian kernel
rbf = @(e,r) exp(-(e*r).^2);

errvec = zeros(length(epvec),1);
dirvec = zeros(length(epvec),1);
likvec = zeros(length(epvec),1);
dlivec = zeros(length(epvec),1);
jvec = zeros(length(epvec),1);
djvec = zeros(length(epvec),1);

k = 1;
for ep=epvec
    GQR = gqr_solve(x, y, ep, alpha);
    errvec(k) = errcompute(gqr_eval(GQR,xx),yy);
    
    b = GQR.coef;
    CbarT = GQR.CbarT;
    lamvec1 = GQR.eig(GQR.Marr(:,1:length(b)));
    lamvec2 = GQR.eig(GQR.Marr(:,length(b)+1:end));
    logdetPsi = sum(log(svd(GQR.stored_psi)));
    logdetPhi = sum(log(svd(GQR.stored_phi1)));
    logdetLam = sum(log(lamvec1));
    logdetK = logdetPsi + logdetPhi + logdetLam;
    bound = b'*(b./lamvec1');
    L2P2P1L1invb = (lamvec2'.^-.5).*(CbarT*b);
    correction = L2P2P1L1invb'*L2P2P1L1invb;
    mahaldist = bound + correction;
    likvec(k) = length(b)*log(mahaldist) + logdetK;
    
    xpf = pick2Dpoints(-1,1,5);
    logdetKtilde = zeros(length(xpf),1);
    for m=1:length(xpf)
        xp = [x;xpf(m,:)];
        GQR2 = gqr_solveprep(0,xp,ep,alpha);
        lamvec = GQR2.eig(GQR2.Marr(:,1:length(xp)));
        Phi1 = GQR2.stored_phi1;
        Phi2 = GQR2.stored_phi2;
        CbarT = GQR2.CbarT;
        logdetPhi = sum(log(svd(Phi1)));
        logdetPsi = sum(log(svd(Phi1+Phi2*CbarT)));
        logdetLam = sum(log(lamvec));
        logdetKtilde(m) = logdetPsi + logdetPhi + logdetLam;
    end
    PF = max(logdetKtilde - logdetK);
    KV_hssvd = log(mahaldist) + log(PF);
    jvec(k) = likvec(k) + KV_hssvd;
    
    K = rbf(ep, DistanceMatrix(x, x));
    Keval = rbf(ep, DistanceMatrix(xx, x));
    warning('off','MATLAB:nearlySingularMatrix')
    c = K\y;
    warning('on','MATLAB:nearlySingularMatrix')
    dirvec(k) = errcompute(Keval*c,yy);
    
    dlivec(k) = length(c)*log(abs(c'*y)) + sum(log(svd(K)));
    
    warning off
    pvals = 1 - sum((Keval/K).*Keval,2);
    warning on
    KV = log(abs(c'*y)) + log(abs(max(pvals)));
    djvec(k) = dlivec(k) + KV;
    
    fprintf('%d, condition of K %e, HS-SVD diff %e\n',k,cond(K),norm(K-GQR.stored_psi*diag(lamvec1)*GQR.stored_phi1'))
    k = k + 1;
end

% h_err = figure;
% loglog(epvec,errvec,'r','linewidth',3)
% hold on
% loglog(epvec,dirvec,'b','linewidth',2)
% hold off
% title('error')
% 
% h_par = figure;
% semilogx(epvec,likvec,'r','linewidth',3)
% hold on
% semilogx(epvec,dlivec,'b','linewidth',2)
% hold off
% title('profile likelihood')
% 
% h_j = figure;
% semilogx(epvec,jvec,'r','linewidth',3)
% hold on
% semilogx(epvec,djvec,'b','linewidth',2)
% hold off
% title('joint metric')

h_full = figure;
[AX,H1,H2] = plotyy(epvec,[dlivec';likvec';djvec';jvec'],epvec,[dirvec';errvec'],'semilogx','loglog');
set(H1,'linewidth',3)
c = get(AX(1),'Children');
set(c(1),'color',[0.8500 0.3250 0.0980])
%set(c(1),'color',[0.9290 0.6940 0.1250])
set(c(1),'linestyle','-')
set(c(2),'color',[0.8500 0.3250 0.0980])
set(c(2),'linestyle',':')
set(c(3),'color',[0 0.4470 0.7410])
set(c(3),'linestyle','-')
set(c(4),'color',[0 0.4470 0.7410])
set(c(4),'linestyle',':')
set(H2,'linewidth',2)
c = get(AX(2),'Children');
set(c(1),'color','k')
set(c(1),'linestyle','-')
set(c(2),'color','k')
set(c(2),'linestyle',':')
set(AX,{'ycolor'},{[.5 0 .5];[0 0 0]}) % Set axis color
set(AX(1),'ylim',[-2500,500])
set(AX(2),'ylim',[1e-14,1e0])
set(AX(2),'ytick',[1e-12,1e-9,1e-6,1e-3])
legend([H1;H2],'MLE direct','MLE HS-SVD','DET direct','DET HS-SVD','Error direct','Error HS-SVD','location','north')
xlabel('\epsilon')
%set(get(AX(1),'xlabel'),'Position',[0.3,-1560,0])
ylabel('Likelihood function')
set(get(AX(1),'ylabel'),'Position',[.005,-1000,17])
set(get(AX(2),'ylabel'),'String','Relative error')
set(get(AX(2),'ylabel'),'Position',[16 3e-07 0])

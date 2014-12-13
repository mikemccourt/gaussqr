% ex17f_gqr.m
% This computes the likelihood for a range of epsilon and fixed input data
% using the HS-SVD for stability in smaller epsilon.
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

NN = 100;
xx = pickpoints(-1,1,NN);
yy = yf(xx);

rbf = @(e,r) exp(-(e*r).^2);
DM_INT = DistanceMatrix(x,x);
DM_EVAL = DistanceMatrix(xx,x);

dirvec = zeros(size(epvec));
gqrvec = zeros(size(epvec));
errvec = zeros(size(epvec));
derrvec = zeros(size(epvec));
boundvec = zeros(size(epvec));
correctionvec = zeros(size(epvec));

k = 1;
h_waitbar = waitbar(0,'Initializing');
for ep=epvec
    % Solve the problem directly
    K = rbf(ep,DM_INT);
    [U,S,V] = svd(K);
    c = (V*diag(0*(diag(S)<1e-14)+(1./diag(S)).*(diag(S)>1e-14))*U'*y);
    K_EVAL = rbf(ep,DM_EVAL);
    derrvec(k) = errcompute(K_EVAL*c,yy);
    
    % Compute the MLE directly
    Mdist = y'*c;
    logdetK = sum(log(diag(S)+eps));
    dirvec(k) = N*log(Mdist) + logdetK;
    
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
    gqrvec(k) = N*log(mahaldist) + logdetK;

    progress = floor(100*k/length(epvec))/100;
    waitbar(progress,h_waitbar,sprintf('Computing, ep=%g',ep))
    k = k + 1;
end

waitbar(100,h_waitbar,sprintf('Plotting'))
h_mle = figure;
[AX,H1,H2] = plotyy(epvec,[gqrvec;dirvec],epvec,[errvec;derrvec],'semilogx','loglog');
set(H1,'linewidth',3)
c = get(AX(1),'Children');
set(c(1),'color',[1 0 0])
set(c(1),'linestyle','--')
set(c(2),'color',[0 0 1])
set(H2,'linewidth',2)
set(H2,'color','k')
set(AX,{'ycolor'},{[.5 0 .5];[0 0 0]}) % Set axis color
set(AX(1),'ylim',[-220,1500])
set(AX(1),'ytick',[0,1000])
set(AX(2),'ylim',[1e-18,1])
set(AX(2),'ytick',[1e-9,1e-5,1e-1])
legend([H1;H2],'Stable MLE','Standard MLE','Stable Error','Standard Error','location','north')
xlabel('\epsilon')
set(get(AX(1),'xlabel'),'Position',[2,-360,0])
ylabel('Likelihood function')
set(get(AX(1),'ylabel'),'Position',[.08,500,17])
set(get(AX(2),'ylabel'),'String','Relative error')
set(get(AX(2),'ylabel'),'Position',[12 3e-014 0])
hold off,pause

h_bound = figure;
gapvec = log10(boundvec./correctionvec);
[AX,H1,H2] = plotyy(epvec,[boundvec;correctionvec],epvec,gapvec,'loglog','semilogx');
set(H1,'linewidth',3)
c = get(AX(1),'Children');
set(c(1),'color',[0 0 0])
set(c(1),'linestyle','--')
set(c(2),'color',[1 0 1])
set(H2,'linewidth',2)
set(H2,'color','b')
set(AX(1),'ylim',[1e-5,1e80])
set(AX(2),'ylim',[-2,3.5])
set(get(AX(2),'xlabel'),'String','\epsilon')
set(get(AX(2),'xlabel'),'Position',[3,-2.23,0])
set(get(H2,'parent'),'ycolor',[0 0 1])
set(AX(1),'ytick',[1,1e20,1e40,1e60])
set(AX(2),'ytick',[-1,0,1,2,3])
legend([H1;H2],'Bound','Correction','Gap','location','southwest')
hold on
minep = min(epvec(gapvec<0));
maxep = max(epvec(gapvec<0));
plot(minep*[1 1],[1e-5,1e80],':','color',[.6,.6,.6],'linewidth',2)
plot(maxep*[1 1],[1e-5,1e80],':','color',[.6,.6,.6],'linewidth',2)
area([minep,maxep],[1e80,1e80],'facecolor',[.9,.9,.9],'basevalue',1e-5)
hold off
uistack(H1(1),'top')
uistack(H1(2),'top') 

close(h_waitbar)
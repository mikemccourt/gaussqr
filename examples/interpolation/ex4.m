% ex4
% This example studies the quality of an interpolatant for different values
% of the ep and alpha parameters required within the Hilbert-Schmidt SVD
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = 200;

epvec = logspace(-3,1,40);
alphavec = logspace(-2,2,39);
N = 45;
NN = 100;

spaceopt = 'cheb';
fopt = 'exp';

rbf = @(e,r) exp(-(e*r).^2);

[yf,fstr] = pickfunc(fopt,1);

[x,spacestr] = pickpoints(-3,3,N,spaceopt);
y = yf(x);
xx = pickpoints(-3,3,NN);
yy = yf(xx);
errvec = zeros(length(epvec),length(alphavec));
errvecd = zeros(size(epvec));
condvec = zeros(length(epvec),length(alphavec));

fprintf(' This could take more than 1 hour to run with M too large \n')

h_waitbar = waitbar(0,'Initializing');
ie = 1;
for ep=epvec
    ia = 1;
    for alpha=alphavec
        GQR = gqr_solve(x,y,ep,alpha);
        yp = gqr_eval(GQR,xx);
        errvec(ie,ia) = errcompute(yp,yy);
        condvec(ie,ia) = strcmp(GQR.warnid,'');
        
        progress = (ia + length(alphavec)*(ie-1))/(length(epvec)*length(alphavec));
        waitbar(progress,h_waitbar,sprintf('ep=%4.3f, alpha=%4.2f',ep,alpha));
        ia = ia + 1;
    end
    
    K = rbf(ep,DistanceMatrix(x,x));
    warning('off','MATLAB:nearlySingularMatrix')
    beta = K\y;
    warning('on','MATLAB:nearlySingularMatrix')
    yp = rbf(ep,DistanceMatrix(xx,x))*beta;
    errvecd(ie) = errcompute(yp,yy);
    
    ie = ie + 1;
end
waitbar(1,h_waitbar,'Generating relevant plots')

[AA,EE] = meshgrid(alphavec,epvec);
errbounded = log10(min(errvec,10^5));

% Plot errors
h_err = figure;

surf(AA,EE,errbounded);
set(gca,'XScale','log')
set(gca,'YScale','log')

ylabel('\epsilon')
xlabel('\alpha')
zlabel('log_{10}(error)')
ptsstr=strcat(', x\in[-3,3],');
title(strcat(fstr,ptsstr,spacestr))
view([.5 -.5 .5])

% Plot condition quality on a contour graph
h_contour = figure;

contour(EE',AA',errbounded');
set(gca,'XScale','log')
set(gca,'YScale','log')

xlabel('\epsilon')
ylabel('\alpha')
title(sprintf('N=%d, x means ill-conditioned',N))
l = colorbar;
set(get(l,'Ylabel'),'String','log_{10}(error)')

badpoints = logical(1-condvec);
text(.95*EE(badpoints),AA(badpoints),'x','FontSize',5)

% Compare results of optimal quality to polynomial results
h_poly = figure;

warning('off','MATLAB:polyfit:RepeatedPoints')
[ppoly,spoly,mupoly] = polyfit(x,y,N-1);
yp = polyval(ppoly,xx,spoly,mupoly);
warning('on','MATLAB:polyfit:RepeatedPoints')
errpoly = norm((yy-yp)./(abs(yy)+eps))/NN;

[minerr,minloc] = min(errbounded,[],2);

[AX,H1,H4] = plotyy(epvec,minerr,epvec,alphavec(minloc),'semilogx','loglog');
hold on
H2 = plot(epvec,log10(errvecd),':xr');
H3 = plot(epvec,log10(errpoly)*ones(size(epvec)),'--k');
xlabel('\epsilon')
set(get(AX(1),'Ylabel'),'String','minimum error')
set(AX(1),'Ylim',[-17,0])
set(AX(1),'YTick',[-15 -10 -5 0])
set(get(AX(2),'Ylabel'),'String','minimizing \alpha')
set(AX(2),'Ylim',[.01,100])
set(AX(2),'Ytick',[.01,1,100])
set(AX(2),'YtickLabel',[.01,1,100])
set(H4,'LineStyle','none')
set(H4,'Marker','o')
legend([H1,H2,H3,H4],'QR','Direct','Polynomial','\alpha','Location','SouthWest')
title(sprintf('f(x)=10e^{-x^2}+x^2, N=%d',N))
hold off

close(h_waitbar)
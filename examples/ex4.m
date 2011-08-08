% Tests different ep and alpha values for a sample problem

epvec = logspace(-3,1,30);
alphavec = logspace(-2,2,29);
N = 30;
NN = 200;

spaceopt = 'cheb';
fopt = 'sinh';

[yf,fstr] = pickfunc(fopt,1);

[x,spacestr] = pickpoints(-3,3,N,spaceopt);
y = yf(x);
xx = pickpoints(-3,3,NN);
yy = yf(xx);
errvec = zeros(length(epvec),length(alphavec));
errvecd = zeros(size(epvec));
condvec = zeros(length(epvec),length(alphavec));

progressbar = 0;
progressincrement = .1;
fprintf(' Progress: ')
ie = 1;
for ep=epvec
    ia = 1;
    for alpha=alphavec
        rbfqrOBJ = rbfqr_solve_alpha(x,y,ep,alpha);
        yp = rbfqr_eval_alpha(rbfqrOBJ,xx);
        errvec(ie,ia) = norm((yy-yp)./(abs(yy)+eps));
        condvec(ie,ia) = strcmp(rbfqrOBJ.warnid,'');
        ia = ia + 1;
    end
    
    K = exp(-ep^2*(repmat(x,1,N)-repmat(x',N,1)).^2);
    warning off MATLAB:nearlySingularMatrix % We know it's bad
    beta = K\y;
    warning on MATLAB:nearlySingularMatrix
    yp = exp(-ep^2*(repmat(x',NN,1)-repmat(xx,1,N)).^2)*beta;
    errvecd(ie) = norm((yy-yp)./(abs(yy)+eps));

    if ie/length(epvec) > progressbar+progressincrement
        progressbar = progressbar + progressincrement;
        fprintf(' %2.0f ',progressbar*100);
    end
    ie = ie + 1;
end
fprintf('\n');

[AA,EE] = meshgrid(alphavec,epvec);
errbounded = log10(min(errvec,10^5));

surf(AA,EE,errbounded);
set(gca,'XScale','log')
set(gca,'YScale','log')

ylabel('\epsilon')
xlabel('\alpha')
zlabel('log_{10}(error)')
ptsstr=strcat(', x\in[-3,3],');
title(strcat(fstr,ptsstr,spacestr))
view([.5 -.5 .5])

figure

contour(EE',AA',errbounded');
set(gca,'XScale','log')
set(gca,'YScale','log')

xlabel('\epsilon')
ylabel('\alpha')
zlabel('log_{10}(error)')
title('Log error contour plot, x means ill-conditioned')
colorbar

badpoints = find(1-condvec);
%hold on
%for k=badpoints % .95 centers the x
    text(.95*EE(badpoints),AA(badpoints),'x','FontSize',5)
%end
%hold off

figure

warning off MATLAB:polyfit:RepeatedPointsOrRescale % We know it's bad
yp = polyval(polyfit(x,y,N-1),xx);
warning on MATLAB:polyfit:RepeatedPointsOrRescale
errpoly = norm((yy-yp)./(abs(yy)+eps));

[minerr,minloc] = min(errbounded,[],2);

[AX,H1,H2] = plotyy(epvec,minerr,epvec,alphavec(minloc),'semilogx','loglog');
hold on
plot(epvec,log10(errvecd),':xr');
plot(epvec,log10(errpoly)*ones(size(epvec)),'--k')
xlabel('\epsilon')
set(get(AX(1),'Ylabel'),'String','minimum error')
set(AX(1),'Ylim',[-13,2])
set(AX(1),'YTick',[-10 -5 0])
set(get(AX(2),'Ylabel'),'String','minimizing \alpha')
set(H2,'LineStyle','none')
set(H2,'Marker','o')
legend('QR','Direct','Polynomial','Location','SouthWest')
title('RBF-QR with optimally chosen \alpha')
hold off
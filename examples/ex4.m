% Tests different ep and alpha values for a sample problem

epvec = logspace(-3,1,40);
alphavec = logspace(-2,2,81);
N = 55;
NN = 200;

spaceopt = 'cheb';
fopt = 'exp';

[yf,fstr] = pickfunc(fopt,1);

[x,spacestr] = pickpoints(-3,3,N,spaceopt);
y = yf(x);
xx = pickpoints(-3,3,NN);
yy = yf(xx);
errvec = zeros(length(epvec),length(alphavec));
errvecd = zeros(size(epvec));
condvec = zeros(length(epvec),length(alphavec));

progressbar = 0;
progressincrement = .03;
fprintf(' Progress: ')
ie = 1;
for ep=epvec
    ia = 1;
    for alpha=alphavec
        rbfqrOBJ = rbfqr_solve_alpha(x,y,ep,alpha);
        yp = rbfqr_eval_alpha(rbfqrOBJ,xx);
        errvec(ie,ia) = norm((yy-yp)./(abs(yy)+eps))/NN;
        condvec(ie,ia) = strcmp(rbfqrOBJ.warnid,'');
        ia = ia + 1;
    end
    
    K = exp(-ep^2*(repmat(x,1,N)-repmat(x',N,1)).^2);
    warning off MATLAB:nearlySingularMatrix % We know it's bad
    beta = K\y;
    warning on MATLAB:nearlySingularMatrix
    yp = exp(-ep^2*(repmat(x',NN,1)-repmat(xx,1,N)).^2)*beta;
    errvecd(ie) = norm((yy-yp)./(abs(yy)+eps))/NN;

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
title(sprintf('N=%d, x means ill-conditioned',N))
l = colorbar;
set(get(l,'Ylabel'),'String','log_{10}(error)')

badpoints = find(1-condvec);
text(.95*EE(badpoints),AA(badpoints),'x','FontSize',5)

figure

warning off MATLAB:polyfit:RepeatedPoints % We know it's bad
[ppoly,spoly,mupoly] = polyfit(x,y,N-1);
yp = polyval(ppoly,xx,spoly,mupoly);
warning on MATLAB:polyfit:RepeatedPoints
errpoly = norm((yy-yp)./(abs(yy)+eps))/NN;

[minerr,minloc] = min(errbounded,[],2);

[AX,H1,H4] = plotyy(epvec,minerr,epvec,alphavec(minloc),'semilogx','loglog');
hold on
H2 = plot(epvec,log10(errvecd),':xr');
H3 = plot(epvec,log10(errpoly)*ones(size(epvec)),'--k');W
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
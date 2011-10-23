% ex2c
% This example compares RBF-Direct to RBF-QR and RBF-QRr for some test
% functions used in optics (see [Jester/Menke/Urban (2011)])

epvecd = logspace(-5,1,20);
epvecr = logspace(-5,1,20);
Nvec = [12,15,18,21;12,15,18,21];
NN = [40;40];

spaceopt = 'even';
fopt = 'KXY';
%fopt = 'KSA1';
%fopt = 'KSA2';
rbf = @(ep,x) exp(-(ep*x).^2);

[yf,fstr] = pickfunc(fopt,2);

%aa = [-1 -1];bb = [1 1];
aa = [-5 -5];bb = [5 5];
%aa = [-20 -20];bb = [20 20];
xx = pick2Dpoints(aa,bb,NN);
yy = yf(xx);

% Use to plot testfunction yf on square grid
[X,Y] = meshgrid(linspace(aa(1),bb(1),40),linspace(aa(2),bb(2),40));
Z = real(yf([X(:),Y(:)]));
Fplot = surfc(X,Y,reshape(Z,40,40));
set(Fplot,'FaceColor','interp','FaceLighting','gouraud','EdgeColor','none')
colormap jet
camlight 
figure

% Use to plot testfunction yf on disk grid
% The following line is not such a good idea to generate points in a disk
% - even though it works - since it reduces NN 
% xx = xx(find(xx(:,1).^2+xx(:,2).^2<=(bb(1)^2)),:);   % points inside disk
% Generate points in a disk using Marco Vianello's wamdisk (weakly
% admissible mesh)
% xx = bb(1)*wamdisk(NN(1,:)-1);
% yy = yf(xx);
% tri = delaunay(xx(:,1),xx(:,2));
% Fplot = trisurf(tri,xx(:,1),xx(:,2),yy);
% set(Fplot,'FaceColor','interp','FaceLighting','gouraud','EdgeColor','none')
% colormap jet
% camlight 
% figure

errvecr = zeros(size(Nvec,2),length(epvecr));
errvecd = zeros(size(Nvec,2),length(epvecd));

status = 'Performing RBF-QRr'
j = 1;
for N=Nvec
    [x,spacestr] = pick2Dpoints(aa,bb,N,spaceopt,ep);
% Use the following lines for approximation on a disk
%    x = bb(1)*wamdisk(N(1)-1);
%    spacestr = 'WAM disk';
    y = yf(x);
    k = 1;
    for ep=epvecr
        rbfqrOBJ = rbfqrr_solve_alpha(x,y,ep);
        yp = rbfqr_eval_alpha(rbfqrOBJ,xx);
        errvecr(j,k) = norm((yy-yp)./(abs(yy)+eps))/sqrt(prod(NN));  % normalized RMS error
        fprintf(' %g ',rbfqrOBJ.alpha)
        %     fprintf(' %d ',k)
        k = k+1;
    end
    fprintf(' %d \n',N(1))
    j = j+1;
end

status = 'Performing RBF-Direct'
j = 1;
for N=Nvec
    [x,spacestr] = pick2Dpoints(aa,bb,N,spaceopt,ep);
% Use the following lines for approximation on a disk
%    x = bb(1)*wamdisk(N(1)-1);
%    spacestr = 'WAM disk';
    y = yf(x);
    k = 1;
    for ep=epvecd
        [x,spacestr] = pick2Dpoints(aa,bb,N,spaceopt,ep);
% Use the following lines for approximation on a disk
%        x = bb(1)*wamdisk(N(1)-1);
%        spacestr = 'WAM disk';
        y = yf(x);
        DM_DATA = DistanceMatrix(x,x);
        IM = rbf(ep,DM_DATA);
        warning off % I know it's bad
        beta = IM\y;
        warning on
        DM_EVAL = DistanceMatrix(xx,x);
        EM = rbf(ep,DM_EVAL);
        yp = EM*beta;
        errvecd(j,k) = norm((yy-yp)./(abs(yy)+eps))/sqrt(prod(NN));  % normalized RMS error
        fprintf(' %d ',k)
        k = k+1;
    end
    fprintf(' %d \n',N(1))
    j = j+1;
end

loglog(epvecd,errvecd,'-.','LineWidth',2)
hold on
loglog(epvecr,errvecr,'LineWidth',3)
hold off
xlabel('\epsilon')
ylabel('normalized RMS error')
ptsstr=strcat(', x\in[',num2str(aa(1)),',',num2str(bb(1)),']^2,');
title(strcat(fstr,ptsstr,spacestr))
legend('N=12^2 (Direct)','N=15^2 (Direct)','N=18^2 (Direct)','N=21^2 (Direct)','N=12^2 (QR)','N=15^2 (QR)','N=18^2 (QR)','N=21^2 (QR)','Location','NorthEast')

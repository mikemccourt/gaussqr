% AppxFit1.m
% This example compares direct RBF interpolation to approximation using the
% HS-SVD eigenfunctions in a regression setting
% The functions of interested come from optical lenses, see
%         [Jester/Menke/Urban (2011)]
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 4;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .3;
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = 100;

epvec = logspace(-2,1,20);
Nvec = [2,4,8,16,32];Nvec = [2,4,8];
NN = 20;

% Choose a kernel, in this case the Gaussian
rbf = @(e,r) exp(-(e*r).^2);

% Choose a function that we want to approximate 
[yf,fstr] = pickfunc('KSA2',2);

% This is just a guess
% In reality, alpha should decrease as ep increases
%   and it should reach a limit as ep->0
alpha = 2;

% Use to plot testfunction yf on disk grid
% The following line is not such a good idea to generate points in a disk
% - even though it works - since it reduces NN 
% xx = xx(find(xx(:,1).^2+xx(:,2).^2<=(bb(1)^2)),:);   % points inside disk
% Generate points in a disk using Marco Vianello's wamdisk (weakly
% admissible mesh)
xx = wamdisk(NN-1);
yy = yf(xx);
tri = delaunay(xx(:,1),xx(:,2));
h_surf = figure;
Fplot = trisurf(tri,xx(:,1),xx(:,2),yy);
set(Fplot,'FaceColor','interp','FaceLighting','gouraud','EdgeColor','none')
colormap jet
camlight

% Run through the GaussQRr regression
h_waitbar = waitbar(0,'Initializing GaussQRr');
j = 1;
errvecr = zeros(size(Nvec,2),length(epvec));
for N=Nvec
    x = wamdisk(N-1);
    y = yf(x);
    k = 1;
    for ep=epvec
        % Remove alpha for an automated orth-based alpha guess
        alpha = max(8*min(1/sqrt(ep),1),2);
        GQR = gqr_solve(x,y,ep,alpha);
        yp = gqr_eval(GQR,xx);
        errvecr(j,k) = errcompute(yp,yy);
        waitbar(((j-1)*length(epvec)+k)/(length(epvec)*size(Nvec,2)),h_waitbar,'Computing GaussQRr')
        k = k+1;
    end
    j = j+1;
end

% Run through the interpolation
waitbar(0,h_waitbar,'Initializing interpolation');
j = 1;
errvecd = zeros(size(Nvec,2),length(epvec));
for N=Nvec
    x = wamdisk(N-1);
    y = yf(x);
    DM = DistanceMatrix(x,x);
    DMeval = DistanceMatrix(xx,x);
    k = 1;
    for ep=epvec
        K = rbf(ep,DM);
        Keval = rbf(ep,DMeval);
        warning('off','MATLAB:nearlySingularMatrix')
        yp = Keval*(K\y);
        warning('on','MATLAB:nearlySingularMatrix')
        errvecd(j,k) = errcompute(yp,yy);
        waitbar(((j-1)*length(epvec)+k)/(length(epvec)*size(Nvec,2)),h_waitbar,'Computing interpolation')
        k = k+1;
    end
    j = j+1;
end

N = 16;
x = wamdisk(N-1);DM = DistanceMatrix(x,x);DMeval = DistanceMatrix(xx,x);
y = yf(x);
alphavec = zeros(size(epvec));
errvec = zeros(size(epvec));
errvecd = zeros(size(epvec));
k = 1;
for ep=epvec
    K = rbf(ep,DM);
    Keval = rbf(ep,DMeval);
    yp = Keval*(K\y);
    errvecd(k) = errcompute(yp,yy);
    [alphavec(k),errvec(k)] = fminbnd(@(alpha)errcompute(gqr_eval(gqr_solve(x,y,ep,alpha),xx),yy),.1,10);
%     errvec(k) = errcompute(gqr_eval(gqr_solve(x,y,ep,alphavec(k)*1.1),xx),yy);
    fprintf('%d\t%g\t%g\t%g\n',k,ep,alphavec(k),errvec(k))
    k = k + 1;
end
figure,loglog(epvec,[errvec;errvecd],'linewidth',3)

waitbar(1,h_waitbar,'Plotting results');

% Plot the observed errors
h_error = figure;
loglog(epvec,errvecd,'-.','LineWidth',2)
hold on
loglog(epvec,errvecr,'LineWidth',3)
hold off
xlabel('\epsilon')
ylabel('RMS relative error')
legendcell = cell(1,2*length(Nvec));
for k=1:length(Nvec)
    legendcell{k} = sprintf('N=%d',Nvec(k)^2+1);
    legendcell{k+length(Nvec)} = sprintf('N=%d',Nvec(k)^2+1);
end
legend(legendcell,'Location','North')

close(h_waitbar)
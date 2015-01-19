% AppxFitOptics.m
% This example studies the use of the low rank eigenfunction approximate
% interpolant as a tool for recovering functions in optics
% The functions of interested come from optical lenses, see
%         [Jester/Menke/Urban (2011)]
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

% Choose a set of data values that we want to test
Nvec = [8,12,16,24,32,40,64];
Neval = 20;

% Choose a kernel, in this case the Gaussian
Mrvec = [.05,.1,.2,.3];
epvec = [.3,1,3];
alpha = 2;

% Choose a function that we want to approximate 
yf = pickfunc('KSA2',2);

% Use to plot testfunction yf on disk grid
% The following line is not such a good idea to generate points in a disk
% - even though it works - since it reduces NN 
% xeval = xeval(find(xeval(:,1).^2+xeval(:,2).^2<=(bb(1)^2)),:);   % points inside disk
% Generate points in a disk using Marco Vianello's wamdisk (weakly
% admissible mesh)
xeval = pick2Dpoints([-1 -1],[1 1],Neval,'wam');
yeval = yf(xeval);
h_wam = figure;
plot(xeval(:,1),xeval(:,2),'.k','linewidth',2)
axis square

tri = delaunay(xeval(:,1),xeval(:,2));
h_surf = figure;
Fplot = trisurf(tri,xeval(:,1),xeval(:,2),yeval);
set(Fplot,'FaceColor','interp','FaceLighting','gouraud','EdgeColor','none')
colormap jet
camlight

% Compute the errors for each of the values of interest
ep = 1;
errmatMr = zeros(length(Mrvec),length(Nvec));
j = 1;
for N=Nvec
    x = pick2Dpoints([-1,-1],[1 1],N,'wam');
    y = yf(x);
    k = 1;
    for Mr=Mrvec
        GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = Mr;
        errmatMr(k,j) = errcompute(gqr_eval(gqr_rsolve(x,y,ep,alpha),xeval),yeval);
        k = k + 1;
    end
    j = j + 1;
end
h_Mrvec = figure;
loglog(Nvec.^2,errmatMr,'linewidth',3)
xlim([min(Nvec.^2),max(Nvec.^2)])
xlabel('number of data points')
ylabel('2-norm error')
legend(cellfun(@(x) ...
    strcat('M=',num2str(x),'N'),num2cell(Mrvec),'UniformOutput',0), ...
    'location','southwest')

% Compute the errors for each of the values of interest
errmatep = zeros(length(Mrvec),length(Nvec));
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .1;
j = 1;
for N=Nvec
    x = pick2Dpoints([-1,-1],[1 1],N,'wam');
    y = yf(x);
    k = 1;
    for ep=epvec
        errmatep(k,j) = errcompute(gqr_eval(gqr_rsolve(x,y,ep,alpha),xeval),yeval);
        k = k + 1;
    end
    j = j + 1;
end
h_epvec = figure;
loglog(Nvec.^2,errmatep,'linewidth',3)
xlim([min(Nvec.^2),max(Nvec.^2)])
xlabel('number of data points')
ylabel('2-norm error')
legend(cellfun(@(x) ...
    strcat('\epsilon=',num2str(x)),num2cell(epvec),'UniformOutput',0), ...
    'location','southwest')
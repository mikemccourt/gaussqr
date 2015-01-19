% AppxFitZonal.m
% This example demonstrates fitting data on a sphere
% Here we will use a zonal kernel rather than a radial kernel
% It may be possible to consider non-zonal kernels on the sphere, but we
% will restrict ourselves to the simpler setting here
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First thing that needs to happen is the Maximal Determinant points need
% to be available on the sphere
% This tests whether these points must be downloaded from the website
% I may move this to rbfsetup at some point
%
% As a disclaimer, these points were taken directly from Rob Wommersley's
% page on Maximal Determinant points
%   http://web.maths.unsw.edu.au/~rsw/Sphere/Extremal/New/extremal1.html
% We thank Dr. Wommersley for his many contributions
gqr_downloaddata('sphereMDpts_data.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Once you have the data, load it into memory
% This will create a cell array sphereMDpts
% These are the points at which we will perform the interpolation
load sphereMDpts_data
Nvec = cellfun(@length,sphereMDpts);

% Choose a test function to interpolate
% First function from Grady Wright at
%   http://math.boisestate.edu/~wright/montestigliano
% Second function from Fasshauer thesis 
%   doi: 10.1016/0377-0427(96)00034-9   [Alfeld, Neamtu, Schumaker]
% yf = @(x) cos(2*(x(:,1)+1/2).^2 + 3*(x(:,2)+1/2).^2 + 5*(x(:,3)-1/sqrt(2)).^2);
yf = @(x) (1 + 10*x(:,1).*x(:,2).*x(:,3) + ...
           x(:,1).^8 + exp(2*x(:,2).^3) + exp(2*x(:,3).^2))/14;

% The original IMQ 1/sqrt(1+(ep*r)^2) is related to the zonal IMQ
% To use the zonal IMQ you can solve the equation
%     gamma/(1-gamma)^2 = ep^2
% with gamma = fzero(@(gamma) gamma./(1-gamma).^2 - ep^2,.5)
% There will always be a gamma in [0,1) for any ep in [0,inf)
% Recall: dp is the dot-product, not the distance
zbf = @(g,dp) 1./sqrt(1+g^2-2*g*dp); % Zonal IMQ
gamma = 0.569; % matches ep = 1.75 for Original IMQ

% Choose various angles and convert to x,y,z
% These are not well distributed, but they are easy
Neval = 200;
xeval = zeros(Neval,3);
[t,testptstr] = pick2Dpoints([-pi -pi/2],[pi pi/2],sqrt(Neval),'halt');
[xeval(:,1),xeval(:,2),xeval(:,3)] = sph2cart(t(:,1),t(:,2),1);
yeval = yf(xeval);

% Study the quality of the interpolation
errvec = zeros(size(sphereMDpts));
for k=1:length(sphereMDpts)
    % Input the data
    x = sphereMDpts{k};
    y = yf(x);
    
    % Create the interpolation and evaluation matrices
    K = zbf(gamma,x*x');
    Ktest = zbf(gamma,xeval*x');
    
    % Perform the interpolation and evaluate the error at the test points
    coef = K\y;
    ytest = Ktest*coef;
    errvec(k) = errcompute(ytest,yeval);
end

% Create some pretty plots of the answer and error
Nplot = 105;
[X,Y,Z] = sphere(Nplot-1);
xplot = [X(:),Y(:),Z(:)];

% This is to plot the true solution
h_true = figure;
yplot = yf(xplot);
C = reshape(yplot,Nplot,Nplot);
surf(X,Y,Z,C,'edgecolor','none')
axis square
view([-1 1 1])
colorbar
colormap gray
title('True solution')

% Now plot the pointwise error in the finest interpolant
% We can reuse the coef from earlier
h_error = figure;
ypplot = zbf(gamma,xplot*x')*coef;
ypdiff = abs(yplot - ypplot);
C = reshape(ypdiff,Nplot,Nplot);
surf(X,Y,Z,C,'edgecolor','none')
view([-1 1 1])
axis square
title(sprintf('Pointwise error, N=%d, gamma=%3.2f',Nvec(end),gamma))
colorbar
colormap gray

% Plot the convergence behavior of this interpolation
h_conv = figure;
loglog(Nvec,errvec,'linewidth',3)
xlabel('N')
ylabel('2-norm error')
title(sprintf('N_{test}=%d%s, gamma=%3.2f',Neval,testptstr,gamma))
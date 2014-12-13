% AppxFit3.m
% This example demonstrates fitting data on a sphere
% Here we will use a zonal kernel rather than a radial kernel
% It may be possible to consider non-zonal kernels on the sphere, but we
% will restrict ourselves to the simpler setting here
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.RANDOM_SEED(0);
dataDir = GAUSSQR_PARAMETERS.DATA_DIRECTORY;
dirslash = GAUSSQR_PARAMETERS.DIRECTORY_SLASH;
GAUSSQR_PARAMETERS.ERROR_STYLE = 4;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First thing that needs to happen is the Maximal Determinant points need
% to be available on the sphere
% This tests whether these points must be downloaded from the website
% I may move this to rbfsetup at some point
%
% As a disclaimer, these points were taken directly from Rob Wommersley's
% page on Maximal Determinant point
%   http://web.maths.unsw.edu.au/~rsw/Sphere/Extremal/New/extremal1.html
% We thank Dr. Wommersley for his many contributions
sphereDataFile = strcat(dataDir,dirslash,'sphereMDpts_data.mat');
if not(exist(sphereDataFile,'file'))
    [filestr,downloadSuccessful] = urlwrite('http://math.iit.edu/~mccomic/gaussqr/data/sphereMDpts_data.mat',sphereDataFile,'Timeout',20);
    if ~downloadSuccessful
        error('Download of data failed or could not be written to %s',filestr)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Once you have the data, load it into memory
% This will create a cell array sphereMDpts
% These are the points at which we will perform the interpolation
load sphereMDpts_data
Nvec = cellfun(@length,sphereMDpts);

% Choose a test function to interpolate
% First function from Grady Wright at
%   http://math.boisestate.edu/~wright/montestigliano/index.html
% Second function from Fasshauer thesis 
%   doi: 10.1016/0377-0427(96)00034-9   [Alfeld, Neamtu, Schumaker]
% yf = @(x) cos(2*(x(:,1)+1/2).^2 + 3*(x(:,2)+1/2).^2 + 5*(x(:,3)-1/sqrt(2)).^2);
yf = @(x) (1 + 10*x(:,1).*x(:,2).*x(:,3) + x(:,1).^8 + exp(2*x(:,2).^3) + exp(2*x(:,3).^2))/14;

% Pick our kernel, here just the inverse multiquadrics
% ep = 2.5 is chosen ad hoc
ep = 2.5;
ZonalDistanceMatrix = @(x,z) sqrt(1-x*z');
zbf = @(e,r) 1./sqrt(1+(e*r).^2);

% Could also choose random or Halton test points around the sphere?
% Choose random angles and convert to x,y,z
% These are not well distributed, but they are easy
NN = 200;
[t,testptstr] = pick2Dpoints([0 0],[pi 2*pi],sqrt(NN),'halt');
xx = [sin(t(:,1)).*cos(t(:,2)),sin(t(:,1)).*sin(t(:,2)),cos(t(:,1))];
yy = yf(xx);

% Study the quality of the interpolation
errvec = zeros(size(sphereMDpts));
for k=1:length(sphereMDpts)
    % Input the data
    x = sphereMDpts{k};
    y = yf(x);
    
    % Create the appropriate distance and interpolation matrices
    DM = ZonalDistanceMatrix(x,x);
    DMtest = ZonalDistanceMatrix(xx,x);
    K = zbf(ep,DM);
    Ktest = zbf(ep,DMtest);
    
    % Perform the interpolation and evaluate the error at the test points
    coef = K\y;
    yp = Ktest*coef;
    errvec(k) = errcompute(yp,yy);
end

% Create some pretty plots of the answer and error
Nplot = 55;
[X,Y,Z] = sphere(Nplot-1);
xplot = [X(:),Y(:),Z(:)];

% This is to plot the true solution
h_true = figure;
yplot = yf(xplot);
C = reshape(yplot,Nplot,Nplot);
surf(X,Y,Z,C,'edgecolor','none')
axis square
view([-1 1 1])
title('True solution')

% Now plot the pointwise error in the finest interpolant
% We can reuse the coef from earlier
h_error = figure;
ypplot = zbf(ep,ZonalDistanceMatrix(xplot,x))*coef;
ypdiff = abs(yplot - ypplot);
C = reshape(ypdiff,Nplot,Nplot);
surf(X,Y,Z,C,'edgecolor','none')
view([-1 1 1])
axis square
colorbar
title(sprintf('Pointwise error, N=%d, ep=%3.2f',Nvec(end),ep))

% Plot the convergence behavior of this interpolation
h_conv = figure;
loglog(Nvec,errvec,'linewidth',3)
xlabel('N')
ylabel('RMS relative error')
title(sprintf('N_{test}=%d%s, ep=%3.2f',NN,testptstr,ep))
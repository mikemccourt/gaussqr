% SurrModelBorehole
% This example creates a surrogate model for the borehole function and
% provides plots of the pairwise errors for each dimension

% Define the function of interest, which is available using the pickfunc
% function available in the GaussQR library
% The _scaled addition allows us to work in [0,1]^8
yf = pickfunc('borehole_scaled');

% If the user has the haltonset file (in Statistics and Machine Learning
% Toolbox), use that to create the Halton points; if not, use the file
% provided in GaussQR
if exist('haltonset','file')
    point_generator = haltonset(8,'Skip',1);
    haltpts8D = @(N) net(point_generator,N);
else
    haltpts8D = @(N) haltonseq(Ntot,8);
end

% Define the kernel for modeling, in its anisotropic form
rbfIM = @(r) 1./(1+r.^2);

% Choose an anisotropic shape parameter vector
% May consider cross-validation at some point
ep = logspace(1,-1,8);

% Create data and test points together
N = 2000;  Neval = 500;
xtot = haltpts8D(N+Neval);

% Identify data and test points and evaluate the function
xeval = xtot(1:Neval,:);  yeval = yf(xeval);
x = xtot(Neval+1:end,:);  y = yf(x);

% Fit the surrogate model and evaluate at the test points
K = rbfIM(DistanceMatrix(x,x,ep));
Keval = rbfIM(DistanceMatrix(xeval,x,ep));
seval = Keval*(K\y);

% Create a vector of pointwise relative errors
errs = (seval - yeval)/norm(yeval)*sqrt(Neval);

% Create the desired contour plots
% To do this, we are required to sample at a regular grid
% We smooth out our interpolant to make the plots easier to see
% We could also use griddata or scatteredInterpolant
muSS = 1e-2;  epSS = 2;
rbfM4 = @(r) (1+r+r.^2/3).*exp(-r);
Npx = 22;  Npy = 21;
xinterp = pick2Dpoints([0 0],[1 1],[Npx Npy]);

% Roll through the dimensions and create the contour plots
% In lieu of the smoothing splines, we could also evaluate the interpolant
% at the baseline values and allow only two dimensions to vary
% This smoothing spline strategy has more of an averaging effect, where
% nearby points with different values have different
% The subplot command starts at the bottom row and walks up
%
% The contour plots and colorbar don't perfectly align, may try to fix it
% at a later time
varnames = {'r_w','r','T_u','H_u','T_\ell','H_\ell','L','K_w'};
h_contour = figure;
hmat = zeros(7);
C = colormap('gray');C = colormap('jet');
for i=1:7
    for j=i:7
        dims = [i,j+1];
        Kdata = rbfM4(DistanceMatrix(xeval(:,dims),xeval(:,dims),epSS));
        Kinterp = rbfM4(DistanceMatrix(xinterp,xeval(:,dims),epSS));
        zinterp = Kinterp*((Kdata + muSS*eye(Neval))\errs);
        
        X = reshape(xinterp(:,1),Npx,Npy);
        Y = reshape(xinterp(:,2),Npx,Npy);
        Z = abs(reshape(zinterp,Npx,Npy));
        subplot(7,7,i+7*(7-j)+(i-1)*7)
        contourf(X,Y,Z,[0,.03,.06,.1],'linewidth',1.5);
        % This colormap command prevents very low errors from looking large
        colormap(.3+C*(.7-.01*min(max(Z(:))/.1,1)));
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'Position',get(gca,'Position').*[1 1 1.25 1.3])
        if i==1
            annotation('textbox',[.06,j/8.17-.02,.1,.1],...
                'string',sprintf('$%s$',varnames{9-j}),...
                'linestyle','none','fontsize',14,'interpreter','latex')
        end
    end
    annotation('textbox',[i/8.6+.03,0,.1,.1],...
        'string',sprintf('$%s$',varnames{i}),...
        'linestyle','none','fontsize',14,'interpreter','latex')
end

% Make the colorbar look good for all the plots
% The ticks and ticklabels changes account for the difference in scaling
% applied to the colormap above
% There might be a better way to do this, but this works
h_color = axes('Position',[.815,.32,0,.62],'Visible','off');
colorbar
allh = get(get(h_color,'parent'),'children');
set(allh(1),'position',get(allh(1),'position').*[1 1 2 1])
set(allh(1),'ticks',[0,1])
set(allh(1),'ticklabels',{'<.03','>.1'})
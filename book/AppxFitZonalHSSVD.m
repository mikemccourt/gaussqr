% AppxFitZonalHSSVD.m
% This script executes the HS-SVD for the spherical harmonics
% Right now, the spherical harmonics evaluation is really terrible
% Again, we STRONGLY recommend that you use Grady Wright's spherical RBF
% introduction http://math.boisestate.edu/~wright/montestigliano if you are
% actually interested in being good at problems on manifolds
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
% We thank Dr. Wommersley for his many contributions to math
gqr_downloaddata('sphereMDpts_data.mat')
load sphereMDpts_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose a test function to interpolate
% First function from Grady Wright at
%   http://math.boisestate.edu/~wright/montestigliano
% Second function from Fasshauer thesis 
%   doi: 10.1016/0377-0427(96)00034-9   [Alfeld, Neamtu, Schumaker]
% yf = @(x) cos(2*(x(:,1)+1/2).^2 + 3*(x(:,2)+1/2).^2 + 5*(x(:,3)-1/sqrt(2)).^2);
yf = @(x) (1 + 10*x(:,1).*x(:,2).*x(:,3) + x(:,1).^8 + exp(2*x(:,2).^3) + exp(2*x(:,3).^2))/14;

% From the data available, use only the N=400 point set for this example
x = sphereMDpts{3};
N = size(x,1);
y = yf(x);
DP = x*x';

% Also choose some points at which to evaluate the error
NN = 2000;
[t,testptstr] = pick2Dpoints([-pi -pi/2],[pi pi/2],sqrt(NN),'halt');
xeval = zeros(NN,3);
[xeval(:,1),xeval(:,2),xeval(:,3)] = sph2cart(t(:,1),t(:,2),1);
yy = yf(xeval);
DPeval = xeval*x';

% The original IMQ 1/sqrt(1+(ep*r)^2) is related to the zonal IMQ
% To use the zonal IMQ you can solve the equation
%     gamma/(1-gamma)^2 = ep^2
% with gamma = fzero(@(gamma) gamma./(1-gamma).^2 - ep^2,.5)
% There will always be a gamma in [0,1) for any ep in [0,inf)
% Recall: dp is the dot-product, not the distance
zbf = @(g,dp) 1./sqrt(1+g^2-2*g*dp);

% Choose a range of shape parameters to test
gammavec = logspace(-2,-.1,30);

% Create the phi matrix, which is full of the first Neig eigenfunctions for
% this kernel.  Although I use the term eigenfunction here, to match the
% notation used for Gaussians, these are the spherical harmonics.  The
% value of Neig would actually need to be a function of epsilon, but if you
% are interested in a more thorough implementation, as opposed to this
% brief example, you should see Grady Wright's work.
% Also, this is done with a loop here for simplicity, but there should
% really be a more methodical way to execute this.
%
% Nlevels is the number of levels of spherical harmonics to consider
% There are l^2-(l-1)^2 = 2*l-1 eigenfunctions in the lth level
% I hate using l as a counter because it looks like 1, but I'm doing it to
% match the notation from various sources
% eigDegree will be needed later for the eigenfunction evaluation
Nlevels = ceil(sqrt(N)) + 3;
Neig = Nlevels^2;
lrange = 0:Nlevels-1;
larr = cell2mat(arrayfun(@(l)l*ones(1,2*l+1),lrange,'UniformOutput',0));
marr = cell2mat(arrayfun(@(l)1:2*l+1,lrange,'UniformOutput',0));
Phi = cell2mat(arrayfun(@(l,m)sphHarm(l,m,x),larr,marr,'UniformOutput',0));
Phieval = cell2mat(arrayfun(@(l,m)sphHarm(l,m,xeval),larr,marr,'UniformOutput',0));
% Break Phi and Phieval into the blocks for the HSSVD
Phi1 = Phi(:,1:N);
Phi2 = Phi(:,N+1:Neig);
Phieval1 = Phieval(:,1:N);
Phieval2 = Phieval(:,N+1:Neig);
Phi2TinvPhi1 = Phi2'/Phi1';

% Perform the interpolation with the range of gamma and record the values
k = 1;
errvec = zeros(size(gammavec));
errvechs = zeros(size(gammavec));
for g=gammavec
    % Perform the standard interpolation
    K = zbf(g,DP);
    Keval = zbf(g,DPeval);
    warning('off','MATLAB:nearlySingularMatrix')
    yeval = Keval*(K\y);
    warning('on','MATLAB:nearlySingularMatrix')
    
    % lamvec is the eigenvalues of the kernel, stored as a vector since
    % storing them as a diagonal matrix seems unnecessary.
    % Also, I realize the 4pi is eliminated automatically, but I am leaving
    % it here for completeness
    % We create LamFull which is the matrix that has the effect of both
    % multiplying on the left by Lam_2 and on the right by inv(Lam_1)
    lamvec = 4*pi*g.^larr./(2*larr+1);
    LamFull = bsxfun(@rdivide,lamvec(N+1:Neig)',lamvec(1:N));
    CbarT = LamFull.*Phi2TinvPhi1;
    Psi = Phi1 + Phi2*CbarT;
    Psieval = Phieval1 + Phieval2*CbarT;
    yevalhs = Psieval*(Psi\y);
    
    errvechs(k) = errcompute(yevalhs,yy);
    errvec(k) = errcompute(yeval,yy);
    k = k + 1;
end

% Plot the results
h = figure;
loglog(gammavec,errvechs,'r','linewidth',3)
hold on
loglog(gammavec,errvec,'--','linewidth',2)
hold off
legend('HS-SVD','Standard Basis','location','northeast')
xlabel('\gamma')
ylabel('2-norm error')
title(sprintf('N=%d points with IMQ kernel',N))
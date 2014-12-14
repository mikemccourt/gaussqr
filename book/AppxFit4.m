% AppxFit4.m
% This script executes the HS-SVD for the spherical harmonics
% Right now, the spherical harmonics evaluation is really terrible
% I would STRONGLY recommend that you consider Grady Wright's spherical RBF
% introduction http://math.boisestate.edu/~wright/montestigliano if you are
% actually interested in being good at problems on manifolds
global GAUSSQR_PARAMETERS
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
gqr_downloaddata('sphereMDpts_data.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose a test function to interpolate
% First function from Grady Wright at
%   http://math.boisestate.edu/~wright/montestigliano
% Second function from Fasshauer thesis 
%   doi: 10.1016/0377-0427(96)00034-9   [Alfeld, Neamtu, Schumaker]
% yf = @(x) cos(2*(x(:,1)+1/2).^2 + 3*(x(:,2)+1/2).^2 + 5*(x(:,3)-1/sqrt(2)).^2);
yf = @(x) (1 + 10*x(:,1).*x(:,2).*x(:,3) + x(:,1).^8 + exp(2*x(:,2).^3) + exp(2*x(:,3).^2))/14;

% Define the distance between points on a sphere
% Note that this is the same as the real Distance Matrix function, it just
% uses a different formulation specific to the sphere
% We write this here to explicitly make this comment
ZonalDistanceMatrix = @(x,z) sqrt(2)*sqrt(1-x*z');

% From the data available, use only the N=400 point set for this example
x = sphereMDpts{3};
N = size(x,1);
y = yf(x);
DM = ZonalDistanceMatrix(x,x);

% Also choose some points at which to evaluate the error
NN = 2000;
[t,testptstr] = pick2Dpoints([-pi -pi/2],[pi pi/2],sqrt(NN),'halt');
xeval = zeros(NN,3);
[xeval(:,1),xeval(:,2),xeval(:,3)] = sph2cart(t(:,1),t(:,2),1);
yy = yf(xeval);
DMeval = ZonalDistanceMatrix(xeval,x);

% The actual RBF in use here is
%     1/sqrt(1+gamma^2 - 2*gamma*x^Tz)
% But I am rewriting it in terms of r, the Euclidean distance
% Note that in this formulation, epsilon is restricted to be within 0 and
% 1, although you can recover the standard IMQ with any such values
rbf = @(e,r) 1./sqrt((1-e)^2+e*r.^2);

% Choose a range of shape parameters to test
epvec = logspace(-2,-.1,30);

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
% There are l^2-(l-1)^2 = 2l-1 eigenfunctions in the lth level
% I hate using l as a counter because it looks like 1, but I'm doing it to
% match the notation from various sources
% eigDegree will be needed later for the eigenfunction evaluation
Nlevels = ceil(sqrt(N)) + 3;
Neig = Nlevels^2;
Phi = zeros(N,Neig);
Phieval = zeros(NN,Neig);
eigDegree = zeros(1,Neig);
k = 1;
for l=0:Nlevels-1
    for m=1:2*l+1
        Phi(:,k) = sphHarm(l,m,x);
        Phieval(:,k) = sphHarm(l,m,xeval);
        eigDegree(k) = l;
        k = k + 1;
    end
end
% Break Phi and Phieval into the blocks for the HSSVD
Phi1 = Phi(:,1:N);
Phi2 = Phi(:,N+1:Neig);
Phieval1 = Phieval(:,1:N);
Phieval2 = Phieval(:,N+1:Neig);
Phi2TinvPhi1 = Phi2'/Phi1';

% Perform the interpolation with the range of ep and record the values
k = 1;
errvec = zeros(size(epvec));
errvechs = zeros(size(epvec));
for ep=epvec
    % Perform the standard interpolation
    K = rbf(ep,DM);
    Keval = rbf(ep,DMeval);
    yeval = Keval*(K\y);
    
    % Eigs are the eigenvalues of the kernel, here stored as a vector since
    % storing them as a diagonal matrix seems unnecessary.
    % Also, I realize the 4pi is eliminated automatically, but I am leaving
    % it here for completeness
    % We create LamFull which is the matrix that has the effect of both
    % multiplying on the left by Lam_2 and on the right by inv(Lam_1)
    Eigs = 4*pi*ep.^eigDegree./(2*eigDegree+1);
    LamFull = repmat(Eigs(N+1:Neig)',1,N)./repmat(Eigs(1:N),Neig-N,1);
    Rbar = LamFull.*Phi2TinvPhi1;
    Psi = Phi1 + Phi2*Rbar;
    Psieval = Phieval1 + Phieval2*Rbar;
    yevalhs = Psieval*(Psi\y);
    
    errvechs(k) = errcompute(yevalhs,yy);
    errvec(k) = errcompute(yeval,yy);
    k = k + 1;
end

% Plot the results
h = figure;
loglog(epvec,errvechs,'r','linewidth',3)
hold on
loglog(epvec,errvec,'--','linewidth',2)
hold off
legend('HS-SVD','Standard Basis','location','north')
xlabel('\epsilon')
ylabel('RMS relative error')
title(sprintf('N=%d points with IMQ kernel',N))
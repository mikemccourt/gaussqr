% ex10.m
% Incorporating the stable HS-SVD basis into MPS
% In honor of Graeme Fairweather's 70th birthday
% This appears in:
%    Using Gaussian eigenfunctions to solve boundary value problems
%    M. McCourt, Advances in Applied Mathematics and Mechanics, 5:569-594, 2013.
%
% The problem we are looking at is 
%     Lap(u) - lambda^2*u = f   -on- interior
%     u = g                     -on- boundary
%   solution: u(x,y) = sin(x^2+y)
%   domain: square between [-1,1]^2
%   fictitious boundary: 1.1 size buffer around real boundary

global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 4;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

% This is the wavenumber (maybe) for the Helmholtz problem
lambda = 3;

% This is the true solution of the problem
fsol = @(x,y) sin(x.^2+y);
% The source associated with the true solution
f = @(x,y) 2*cos(x.^2+y)-4*x.^2.*sin(x.^2+y)-sin(x.^2+y)-lambda^2*fsol(x,y);

% An alternate problem
% fsol = @(x,y) exp(x+y);
% f = @(x,y) (2-lambda^2)*exp(x+y);

% The fundamental solution of the Helmholtz problem
Hfs = @(r) besselk(0,lambda*r)/(2*pi);

% The number of error evalution points in each dimension
NN = 35;

% The length of the GaussQRr regression
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .5;
% The global scale parameter for the HS-SVD basis
alpha = 1;
% The Gaussian shape parameter
ep = 1e-5;

bvec = 10:5:80;
errMPS = zeros(size(bvec));
errMPSimp = zeros(size(bvec));
errGQR = zeros(size(bvec));
bNvec = zeros(size(bvec));

m = 1;
for bN=bvec
    % We need both the collocation and the source points for MFS
    ptsMFScoll = [[linspace(-1,1,bN)';linspace(-1,1,bN)';-ones(bN,1);ones(bN,1)],...
        [-ones(bN,1);ones(bN,1);linspace(-1,1,bN)';linspace(-1,1,bN)']];
    ptsMFScoll = unique(1e-8*ceil(1e8*ptsMFScoll),'rows');
    bNvec(m) = size(ptsMFScoll,1);
    % Use a circle around the boundary (if preferred)
%     ptsMFSsource = 1.5*[cos(linspace(-pi,pi,size(ptsMFScoll,1)))',sin(linspace(-pi,pi,size(ptsMFScoll,1)))'];
    % Use points spaced close to the boundary
    ptsMFSsource = 1/lambda^2*ptsMFScoll.*(abs(ptsMFScoll)==1) + ptsMFScoll;
    
    % We also need to set up the interior solution by GaussQR
    ptsGQR = pick2Dpoints([-1 -1],[1 1],floor(sqrt(bNvec(m))),'halton');
    ptsGQR = 1e-8*ceil(1e8*ptsGQR);
    GQRonBDY = find(any(abs(ptsGQR)==1,2)); % find points on the boundary
    % Dump points on the boundary
    ptsGQR = ptsGQR(setdiff(1:size(ptsGQR,1),GQRonBDY),:);
    % Restrict GQR to only using the same # of points at MFS or less
    ptsGQR = ptsGQR(setdiff(1:min(size(ptsGQR,1),bNvec(m)),GQRonBDY),:);
    
    % Find some sample points to evaluate the error at
    ptsEVAL = pick2Dpoints([-1 -1],[1 1],NN);
    % Evaluate the true solution at the sample points
    usol = fsol(ptsEVAL(:,1),ptsEVAL(:,2));
    
    % Solve the particular solution problem indirectly
    GQR = gqr_solveprep(1,ptsGQR,ep,alpha);
    phiMat = gqr_phi(GQR,ptsGQR);
    phiMat2d = gqr_phi(GQR,ptsGQR,[2,0])+gqr_phi(GQR,ptsGQR,[0,2]);
    A = phiMat2d - lambda^2*phiMat;
    rhs = f(ptsGQR(:,1),ptsGQR(:,2));
    
    GQR.coef = A\rhs;
    
    % Store the size of the Marr for use below, to make it fair
    M_MPS = size(GQR.Marr,2);
    
    % Consider the problem with just GaussQR for comparison
    % This line allows for more boundary points as N increases
    ptsBDY = pick2Dpoints([-1,-1],[1 1],sqrt(size(ptsGQR,1)));
    ptsBDY = ptsBDY(any(abs(ptsBDY)==1,2),:);
    ptsFULL = [ptsGQR;ptsBDY];
    
    % Form the solution using only the HS-SVD solution
    GQRfull = gqr_solveprep(1,ptsFULL,ep,alpha,M_MPS);
    phiMat = gqr_phi(GQRfull,ptsGQR);
    phiMatBC = gqr_phi(GQRfull,ptsBDY);
    phiMat2d = gqr_phi(GQRfull,ptsGQR,[2,0])+gqr_phi(GQRfull,ptsGQR,[0,2]);
    A = [phiMat2d - lambda^2*phiMat;phiMatBC];
    rhs = [f(ptsGQR(:,1),ptsGQR(:,2));fsol(ptsBDY(:,1),ptsBDY(:,2))];
    
    GQRfull.coef = A\rhs;
    
    % Now enforce the boundary with MFS
    uPonBDY = gqr_eval(GQR,ptsMFScoll);
    rhs = fsol(ptsMFScoll(:,1),ptsMFScoll(:,2)) - uPonBDY;
    DM_coll = DistanceMatrix(ptsMFScoll,ptsMFSsource);
    A_coll = Hfs(DM_coll);
    coefMFS = A_coll\rhs;
    
    % Also consider improved MFS with full GaussQR
    uPonBDY = gqr_eval(GQRfull,ptsMFScoll);
    rhs = fsol(ptsMFScoll(:,1),ptsMFScoll(:,2)) - uPonBDY;
    DM_coll = DistanceMatrix(ptsMFScoll,ptsMFSsource);
    A_coll = Hfs(DM_coll);
    coefMFSimp = A_coll\rhs;
    
    % Now evaluate the full solution at the error points
    uP_eval = gqr_eval(GQR,ptsEVAL);
    DM_eval = DistanceMatrix(ptsEVAL,ptsMFSsource);
    A_eval = Hfs(DM_eval);
    uF_eval = A_eval*coefMFS;
    uPi_eval = gqr_eval(GQRfull,ptsEVAL);
    uFi_eval = A_eval*coefMFSimp;
    errMPS(m) = errcompute(uF_eval+uP_eval,usol);
    errGQR(m) = errcompute(gqr_eval(GQRfull,ptsEVAL),usol);
    errMPSimp(m) = errcompute(uFi_eval+uPi_eval,usol);
    m = m+1;
end

% Plot of error with MPSimp
loglog(bNvec,[errMPS;errGQR;errMPSimp],'linewidth',3)
ylabel('RMS error')
xlabel('Collocation points')
xlim([min(bNvec),max(bNvec)])
legend('MPS','GaussQR','MPS+GaussQR','location','southwest')
set(gca,'xtick',[36,50,80,140,250])
set(gca,'ytick',[1e-12,1e-9,1e-6,1e-3])
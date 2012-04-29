% ex5d.m
% Incorporating GaussQR into MPS
% In honor of Graeme Fairweather's 70th birthday
% The first problem we'll be looking at is 
%     Lap(u) - lambda^2*u = f   -on- interior
%     u = g                     -on- boundary
%   solution: u(x,y) = exp(x+y)
%   source:   f(x,y) = (2-lambda)exp(x+y)
%   domain: square between [-1,1]^2
%   fictitious boundary: 1.1 times real boundary

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
GAUSSQR_PARAMETERS.ERROR_STYLE = 3;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

% This is the wavenumber (maybe) for the Helmholtz problem
lambda = 3;

% This is the true solution of the problem
fsol = @(x,y) sin(x.^2+y);
% The source associated with the true solution
f = @(x,y) 2*cos(x.^2+y)-4*x.^2.*sin(x.^2+y)-sin(x.^2+y)-lambda^2*fsol(x,y);

fsol = @(x,y) exp(x+y);
f = @(x,y) (2-lambda^2)*exp(x+y);
% The fundamental solution of the Helmholtz problem
Hfs = @(r) besselk(0,lambda*r)/(2*pi);

% The number of error evalution points in each dimension
NN = 35;

% The length of the GaussQRr regression
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .6;
% The global scale parameter for GaussQR
alpha = 1;
% The shape parameter commonly associated with RBFs
ep = 1e-9;

bvec = 10:5:70;
errMPS = [];
errMPSimp = [];
errGQR = [];
bNvec = [];

m = 1;
for bN=bvec
    % We need both the collocation and the source points for MFS
    ptsMFScoll = [[linspace(-1,1,bN)';linspace(-1,1,bN)';-ones(bN,1);ones(bN,1)],...
        [-ones(bN,1);ones(bN,1);linspace(-1,1,bN)';linspace(-1,1,bN)']];
    ptsMFScoll = unique(1e-8*ceil(1e8*ptsMFScoll),'rows');
    bNvec(m) = size(ptsMFScoll,1);
    % Use a circle around the boundary
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
    % ptsGQR = ptsGQR(setdiff(1:min(size(ptsGQR,1),bNvec(m)),GQRonBDY),:);
    
    % Find some sample points to evaluate the error at
    ptsEVAL = pick2Dpoints([-1 -1],[1 1],NN);
    % Evaluate the true solution at the sample points
    usol = fsol(ptsEVAL(:,1),ptsEVAL(:,2));
    
    % Solve the particular solution problem indirectly
    [ep,alpha,Marr] = gqr_solveprep(1,ptsGQR,ep,alpha);
    phiMat = gqr_phi(Marr,ptsGQR,ep,alpha);
    phiMat2d = gqr_phi(Marr,ptsGQR,ep,alpha,[2,0])+gqr_phi(Marr,ptsGQR,ep,alpha,[0,2]);
    A = phiMat2d - lambda^2*phiMat;
    rhs = f(ptsGQR(:,1),ptsGQR(:,2));
    coef = A\rhs;
    
    % Fill the GQR object with all the values it needs
    GQR.reg = 1;
    GQR.Marr = Marr;
    GQR.alpha = alpha;
    GQR.ep = ep;
    GQR.N = size(ptsGQR,1);
    GQR.coef = coef;
    
    % Consider the problem with just GaussQR for comparison
    % only a fixed number of boundary points are used here
    ptsBDY = pick2Dpoints([-1,-1],[1 1],6);
    % This line allows for more boundary points as N increases
%     ptsBDY = pick2Dpoints([-1,-1],[1 1],sqrt(GQR.N));
    ptsBDY = ptsBDY(find(any(abs(ptsBDY)==1,2)),:);
    ptsFULL = [ptsGQR;ptsBDY];
    [ep,alpha,Marr] = gqr_solveprep(1,ptsFULL,ep,alpha);
    phiMat = gqr_phi(Marr,ptsGQR,ep,alpha);
    phiMatBC = gqr_phi(Marr,ptsBDY,ep,alpha);
    phiMat2d = gqr_phi(Marr,ptsGQR,ep,alpha,[2,0])+gqr_phi(Marr,ptsGQR,ep,alpha,[0,2]);
    A = [phiMat2d - lambda^2*phiMat;phiMatBC];
    rhs = [f(ptsGQR(:,1),ptsGQR(:,2));fsol(ptsBDY(:,1),ptsBDY(:,2))];
    coef = A\rhs;
    
    % Fill the GQR object with all the values it needs
    GQRfull.reg = 1;
    GQRfull.Marr = Marr;
    GQRfull.alpha = alpha;
    GQRfull.ep = ep;
    GQRfull.N = size(ptsFULL,1);
    GQRfull.coef = coef;
    
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
    errGQR(m) = errcompute(rbfqr_eval(GQRfull,ptsEVAL),usol);
    errMPSimp(m) = errcompute(uFi_eval+uPi_eval,usol);
    m = m+1;
end

% Plot of error with MFS
loglog(bNvec,[errMPS;errGQR;errMPSimp])
ylabel('error')
xlabel('Collocation points')
title('Solution via MPS')
xlim([min(bNvec),max(bNvec)])
legend('MPS','GaussQR','MPSimp')
%set(gca,'xtick',bNvec)

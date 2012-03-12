% ex5d.m
% Comparing GaussQR to MFS for a couple problems
% In honor of Graeme Fairweather's 70th birthday
% The first problem we'll be looking at is 
%     Lap(u) = 0
%     Dirichlet BC
%   solution: u(x,y) = exp(x)cos(y)
%   domain: rectangle between (0,0)->(1,pi/2)
%   fictitious boundary: Circle, radius 2, center (.5,pi/4)

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
GAUSSQR_PARAMETERS.ERROR_STYLE = 3; % Use absolute error

fsol = @(x,y) exp(x).*cos(y);
Lfs = @(r) -1/(2*pi)*log(r);
NN = 20;
bvec = 6:20;
errMFS = [];

m = 1;
for bN=bvec
    % We need both the collocation and the source points for MFS
    ptsMFScoll = [[linspace(0,1,bN)';linspace(0,1,bN)';zeros(bN,1);ones(bN,1)],...
        [zeros(bN,1);pi/2*ones(bN,1);linspace(0,pi/2,bN)';linspace(0,pi/2,bN)']];
    ptsMFScoll = unique(ptsMFScoll,'rows');
    ptsMFSsource = 2*[cos(linspace(-pi,pi,size(ptsMFScoll,1)))',sin(linspace(-pi,pi,size(ptsMFScoll,1)))'] + ones(size(ptsMFScoll,1),1)*[.5,pi/4];
    % Find some sample points to evaluate the error at
    ptsEVAL = pick2Dpoints([0 0],[1 pi/2],NN);
    % Evaluate the true solution at the sample points
    usol = fsol(ptsEVAL(:,1),ptsEVAL(:,2));

    % Solve the system first with MFS
    rhs = fsol(ptsMFScoll(:,1),ptsMFScoll(:,2));
    DM_coll = DistanceMatrix(ptsMFScoll,ptsMFSsource);
    A_coll = Lfs(DM_coll);
    coefMFS = A_coll\rhs;
    DM_eval = DistanceMatrix(ptsEVAL,ptsMFSsource);
    A_eval = Lfs(DM_eval);
    uMFS = A_eval*coefMFS;
    errMFS(m) = errcompute(uMFS,usol);
    m = m + 1;
end

% Contour plot of error with MFS
subplot(1,2,1)
loglog(4*(bvec-1),errMFS)
ylabel('error')
xlabel('Collocation points')
title('Solution via MFS')
set(gca,'xtick',4*(bvec-1))

% Solve the system with GaussQRr with full regression
fsol = @(x,y) exp(.5*x+.5).*cos(pi/4*(y+1));
f = @(x,y) (1/4-pi^2/16)*fsol(x,y);

ptsEVAL = pick2Dpoints([-1 -1],[1 1],NN);
usol = fsol(ptsEVAL(:,1),ptsEVAL(:,2));

GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .8;
alpha = 2;
ep = 1e-9;
Marr = rbfformMarr([0;0],[],length(ptsRBF));

m = 1;
errvecR2D = [];
Nvec = [];
for N=5:10
    x = unique([pick2Dpoints([-1 -1],[1 1],N,'cheb');pick2Dpoints([-1 -1],[1 1],N,'halton')],'rows');
    b = find(abs(x(:,1))==1 | abs(x(:,2))==1);
    bi = setdiff(1:size(x,1),b)';
    Nvec(m) = size(x,1);
    
    [ep,alpha,Marr] = rbfsolveprep(1,x,ep,alpha);
    phiMat = rbfphi(Marr,x(b,:),ep,alpha);
    phiMat2d = rbfphi(Marr,x(bi,:),ep,alpha,[2,0])+rbfphi(Marr,x(bi,:),ep,alpha,[0,2]);
    A = [phiMat2d;phiMat];
    rhs = zeros(length(x),1);
    rhs(1:length(bi)) = f(x(bi,1),x(bi,2));
    rhs(1+length(bi):end) = fsol(x(b,1),x(b,2));
    coef = A\rhs;
    
    GQR.reg = 1;
    GQR.Marr = Marr;
    GQR.alpha = alpha;
    GQR.ep = ep;
    GQR.N = length(ptsRBF);
    GQR.coef = coef;
    
    uGQR = rbfqr_eval(GQR,ptsEVAL);
    errvecR2D(m) = errcompute(uGQR,usol);
    m = m+1;
end

subplot(1,2,2)
loglog(Nvec,errvecR2D)
ylabel('error')
xlabel('Collocation points')
title('Solution via GaussQR')
xlim([min(Nvec),max(Nvec)])
ylim([1e-16 1e-6])
set(gca,'xtick',Nvec)
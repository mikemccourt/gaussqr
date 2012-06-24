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
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

fsol = @(x,y) exp(x).*cos(y);
Lfs = @(r) -1/(2*pi)*log(r);
NN = 25;
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .8;

bvec = 4:15;
errMFS = [];

m = 1;
for bN=bvec
    % We need both the collocation and the source points for MFS
    ptsMFScoll = [[linspace(0,1,bN)';linspace(0,1,bN)';zeros(bN,1);ones(bN,1)],...
        [zeros(bN,1);pi/2*ones(bN,1);linspace(0,pi/2,bN)';linspace(0,pi/2,bN)']];
    ptsMFScoll = unique(1e-8*ceil(1e8*ptsMFScoll),'rows');
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
NvecMFS = 4*(bvec-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the system with GaussQRr with partial regression
fsol = @(x,y) exp(.5*x+.5).*cos(pi/4*(y+1));
f = @(x,y) (1/4-pi^2/16)*fsol(x,y);

ptsEVAL = pick2Dpoints([-1 -1],[1 1],NN);
usol = fsol(ptsEVAL(:,1),ptsEVAL(:,2));

alpha = 1;
ep = 1e-8;

m = 1;
errGQR = [];
Nvec = [];
for N=3:10
    % Choose some possible collocation points and eliminate duplicates
    x = [pick2Dpoints([-1 -1],[1 1],N,'cheb');pick2Dpoints([-1 -1],[1 1],N,'halton')];
    x = unique(1e-8*ceil(1e8*x),'rows');
    b = find(abs(x(:,1))==1 | abs(x(:,2))==1);
    bi = setdiff(1:size(x,1),b)';
    Nvec(m) = size(x,1);
    
    [ep,alpha,Marr] = gqr_solveprep(1,x,ep,alpha);
    phiMat = gqr_phi(Marr,x(b,:),ep,alpha);
    phiMat2d = gqr_phi(Marr,x(bi,:),ep,alpha,[2,0])+gqr_phi(Marr,x(bi,:),ep,alpha,[0,2]);
    A = [phiMat2d;phiMat];
    rhs = zeros(length(x),1);
    rhs(1:length(bi)) = f(x(bi,1),x(bi,2));
    rhs(1+length(bi):end) = fsol(x(b,1),x(b,2));
    coef = A\rhs;
    
    GQR.reg = 1;
    GQR.Marr = Marr;
    GQR.alpha = alpha;
    GQR.ep = ep;
    GQR.N = length(x);
    GQR.coef = coef;
    
    uGQR = gqr_eval(GQR,ptsEVAL);
    errGQR(m) = errcompute(uGQR,usol);
    m = m+1;
end
NvecGQR = Nvec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider a 4th order FD discretization
fsol = @(x) exp(.5*x(:,1)+.5).*cos(pi/4*(x(:,2)+1));
f = @(x) (1/4-pi^2/16)*fsol(x);
stencil9 = [1/6 2/3 1/6  2/3 -10/3 2/3  1/6 2/3 1/6];

Nvec = 4:20;
N2vec = [];
errFD = [];
k = 1;
for N=Nvec
    x = pick2Dpoints([-1,-1],[1,1],N);
    dx = 2/(N-1);
    usol = fsol(x);
    
    A = sparse(N,N);
    b = find(abs(x(:,1))==1 | abs(x(:,2))==1)';
    bi = setdiff(1:size(x,1),b);
    Astenoff9 = [-N-1 -N -N+1  -1 0 1  N-1 N N+1];
    for ix=bi
        A(ix,ix+Astenoff9) = stencil9;
    end
    A(b,b) = speye(length(b)); % Handle dirichlet BC
    
    rhs = zeros(N,1);
    rhs(bi) = dx^2*f(x(bi,:));
    rhs(b) = fsol(x(b,:));
    
    uFD = A\rhs;
    
    N2vec(k) = size(x,1);
    errFD(k) = errcompute(uFD,usol);
    k = k+1;
end
NvecFD = N2vec;

loglog(NvecMFS,errMFS,'b','linewidth',3),hold on
loglog(NvecGQR,errGQR,'k','linewidth',3)
loglog(NvecFD,errFD,'r','linewidth',3),hold off
xlabel('Number of solution points')
ylabel('Absolute inf-norm error')
legend('MFS','GaussQRr','4th order FD','location','east')
xlim([16,200])
ylim([1e-15 1e-1])
set(gca,'xtick',[16,25,50,100,200])
% ex13b_gqr.m
% Comparing GaussQR to MFS for a couple problems
% In honor of Graeme Fairweather's 70th birthday
% The first problem we'll be looking at is 
%     Lap(u) = 0
%   solution: u(x,y) = exp(x)cos(y)
%   domain: L-shaped domain between (0,0)->(1,pi/2)
%   fictitious boundary: Circle, radius 2, center (.5,pi/4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The boundary conditions are chosen to be more difficult
%   On x=0, we use the Neumann BC
%   Everywhere else, we use the Dirichlet BC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;
fsol = @(x,y) exp(x).*cos(y);
fdysol = @(x,y) -exp(x).*sin(y);
Lfs = @(r) -1/(2*pi)*log(r);
Ldyfs = @(x,y) -1./(2*pi*DistanceMatrix(x,y).^2).*DifferenceMatrix(x(:,2),y(:,2));

% First, we're gonna look at the Finite Difference solution
bvec = 9:20;
bNvec = [];
errFD = [];

m = 1;
for bN=bvec
    % Find the computational grid
    L = numgrid('L',bN)';
    N = max(max(L));
    bNvec(m) = N
    A = delsq(L);
    BCpts = find(full(sum(delsq(L)~=0,2))<5);
    BCNpts = L(end-1,2:end-1);
    BCDpts = setdiff(BCpts,BCNpts);
    INTpts = setdiff(1:N,BCpts);
    
    % Fill in the solution matrix
    dx = pi/2/(N-1);
    dy = 1/(N-1);
    for k=INTpts
        xk = find(L==k);
        Lc = ceil(xk/bN);
        Lr = mod(xk,bN);
        A(xk*ones(1,5),[L(Lr-1,Lc),L(Lr,Lc-1),L(Lr,Lc),L(Lr,Lc+1),L(Lr+1,Lc)]) = ...
             [-1/dy^2,-1/dx^2,2/dy^2+2/dx^2,-1/dx^2,-1/dy^2];
    end
end

pause

bvec = 5:10:55;
bNvec = [];
errMFS = [];

m = 1;
for bN=bvec
    % Choose the collocation points
    x = [pick2Dpoints([-1 -1],[1 1],bN);pick2Dpoints([0 0],[1 1],ceil(bN/2))];
    bx = find( x(:,1)==-1 | x(:,2)==-1 | (x(:,1)==1 & x(:,2)<=0) | (x(:,2)==1 & x(:,1)<=0) | (x(:,1)>=0 & x(:,2)==0) | (x(:,2)>=0 & x(:,1)==0) );
    x = unique(1e-8*ceil(1e8*x(bx,:)),'rows');
    ptsMFScoll = x*diag([.5 pi/4]) + ones(size(x,1),1)*[.5 pi/4];
    bNvec(m) = size(ptsMFScoll,1);

    % Find the points which are Neumann BC instead of Dirichlet
    b = find( ptsMFScoll(:,2)==0 );
    
    % Choose the source points, just the standard circle
    ptsMFSsource = 2*[cos(linspace(-pi,pi,bNvec(m)))',sin(linspace(-pi,pi,bNvec(m)))'] + ones(bNvec(m),1)*[.5,pi/4];

    % Find some sample points at which to evaluate the error
    x = pick2Dpoints([-1 -1],[1 1],NN);
    bx = find( x(:,1)<=0 | x(:,2)<=0 );
    ptsEVAL = x(bx,:)*diag([.5 pi/4]) + ones(size(bx,1),1)*[.5 pi/4];
    % Evaluate the true solution at the sample points
    usol = fsol(ptsEVAL(:,1),ptsEVAL(:,2));

    % Solve the system with MFS
    DM_coll = DistanceMatrix(ptsMFScoll,ptsMFSsource);
    % Handle the Dirichlet BC
    A_coll = Lfs(DM_coll);
    rhs = fsol(ptsMFScoll(:,1),ptsMFScoll(:,2));
    % Handle the Neumann BC
    A_coll(b,:) = Ldyfs(ptsMFScoll(b,:),ptsMFSsource);
    rhs(b) = fdysol(ptsMFScoll(b,1),ptsMFScoll(b,2));

    coefMFS = A_coll\rhs;
    DM_eval = DistanceMatrix(ptsEVAL,ptsMFSsource);
    A_eval = Lfs(DM_eval);
    uMFS = A_eval*coefMFS;
    errMFS(m) = errcompute(uMFS,usol);
    m = m + 1;
end


% GaussQR solution
fsol = @(x,y) exp(.5*x+.5).*cos(pi/4*(y+1));
fdysol = @(x,y) -pi/4*exp(.5*x+.5).*sin(pi/4*(y+1));
f = @(x,y) (1/4-pi^2/16)*fsol(x,y);

ptsEVAL = pick2Dpoints([-1 -1],[1 1],NN);
ptsEVAL = ptsEVAL(find( ptsEVAL(:,1)<=0 | ptsEVAL(:,2)<=0 ),:);
usol = fsol(ptsEVAL(:,1),ptsEVAL(:,2));

alpha = 1;
ep = 1e-9;

m = 1;
errvecR2D = [];
Nvec = [];
for N=5:12
    % Determine the collocation points
    x = [pick2Dpoints([-1 -1],[1 1],N,'cheb');pick2Dpoints([-1 -1],[1 1],N,'halton');pick2Dpoints([0 0],[1 1],ceil(N/2))];
    x = unique(1e-8*ceil(1e8*x),'rows');
    bx = find( x(:,1)<=0 | x(:,2)<=0 );
    x = x(bx,:);
    Nvec(m) = size(x,1);

    % Find out which points are on the Dirichlet boundary
    bD = find( abs(x(:,1))==1 | x(:,2)==1 | (x(:,1)==0 & x(:,2)>=0) | (x(:,1)>=0 & x(:,2)==0));
    % Find out which points are on the Neumann boundary
    bN = find( x(:,2)==-1 &  abs(x(:,1))~=1 );
    % Find out which points are on the interior
    bi = setdiff(1:Nvec(m),[bN;bD])';

    % Get the preliminary RBF-QR stuff (namely Marr)
    [ep,alpha,Marr] = gqr_solveprep(1,x,ep,alpha);
    % Build the collocation matrix
    phiMat2d = gqr_phi(Marr,x(bi,:),ep,alpha,[2,0])+gqr_phi(Marr,x(bi,:),ep,alpha,[0,2]);
    phiMat = gqr_phi(Marr,x(bD,:),ep,alpha);
    phiMat1dy = gqr_phi(Marr,x(bN,:),ep,alpha,[0,1]);
    A = [phiMat2d;phiMat;phiMat1dy];
    % Build the RHS
    rhs_interior = f(x(bi,1),x(bi,2));
    rhs_dirichlet = fsol(x(bD,1),x(bD,2));
    rhs_neumann = fdysol(x(bN,1),x(bN,2));
    rhs = [rhs_interior;rhs_dirichlet;rhs_neumann];
    coef = A\rhs;

    GQR.reg = 1;
    GQR.Marr = Marr;
    GQR.alpha = alpha;
    GQR.ep = ep;
    GQR.N = Nvec(m);
    GQR.coef = coef;

    uGQR = gqr_eval(GQR,ptsEVAL);
    errvecR2D(m) = errcompute(uGQR,usol);
    m = m + 1;
end


%Finite Difference solution
fsol = @(x,y) exp(.5*x+.5).*cos(pi/4*(y+1));
fdysol = @(x,y) -pi/4*exp(.5*x+.5).*sin(pi/4*(y+1));
f = @(x,y) (1/4-pi^2/16)*fsol(x,y);

ptsEVAL = pick2Dpoints([-1 -1],[1 1],NN);
ptsEVAL = ptsEVAL(find( ptsEVAL(:,1)<=0 | ptsEVAL(:,2)<=0 ),:);
usol = fsol(ptsEVAL(:,1),ptsEVAL(:,2));

alpha = 1;
ep = 1e-9;

m = 1;
errvecR2D = [];
Nvec = [];
for N=5:12
    % Determine the collocation points
    x = [pick2Dpoints([-1 -1],[1 1],N,'cheb');pick2Dpoints([-1 -1],[1 1],N,'halton');pick2Dpoints([0 0],[1 1],ceil(N/2))];
    x = unique(1e-8*ceil(1e8*x),'rows');
    bx = find( x(:,1)<=0 | x(:,2)<=0 );
    x = x(bx,:);
    Nvec(m) = size(x,1);

    % Find out which points are on the Dirichlet boundary
    bD = find( abs(x(:,1))==1 | x(:,2)==1 | (x(:,1)==0 & x(:,2)>=0) | (x(:,1)>=0 & x(:,2)==0));
    % Find out which points are on the Neumann boundary
    bN = find( x(:,2)==-1 &  abs(x(:,1))~=1 );
    % Find out which points are on the interior
    bi = setdiff(1:Nvec(m),[bN;bD])';

    % Get the preliminary RBF-QR stuff (namely Marr)
    [ep,alpha,Marr] = gqr_solveprep(1,x,ep,alpha);
    % Build the collocation matrix
    phiMat2d = gqr_phi(Marr,x(bi,:),ep,alpha,[2,0])+gqr_phi(Marr,x(bi,:),ep,alpha,[0,2]);
    phiMat = gqr_phi(Marr,x(bD,:),ep,alpha);
    phiMat1dy = gqr_phi(Marr,x(bN,:),ep,alpha,[0,1]);
    A = [phiMat2d;phiMat;phiMat1dy];
    % Build the RHS
    rhs_interior = f(x(bi,1),x(bi,2));
    rhs_dirichlet = fsol(x(bD,1),x(bD,2));
    rhs_neumann = fdysol(x(bN,1),x(bN,2));
    rhs = [rhs_interior;rhs_dirichlet;rhs_neumann];
    coef = A\rhs;

    GQR.reg = 1;
    GQR.Marr = Marr;
    GQR.alpha = alpha;
    GQR.ep = ep;
    GQR.N = Nvec(m);
    GQR.coef = coef;

    uGQR = gqr_eval(GQR,ptsEVAL);
    errvecR2D(m) = errcompute(uGQR,usol);
    m = m + 1;
end


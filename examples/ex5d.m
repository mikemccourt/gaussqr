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
GAUSSQR_PARAMETERS.ERROR_STYLE = 3; % Use partial relative error

fsol = @(x,y) exp(x).*cos(y);
Lfs = @(r) -1/(2*pi)*log(r);
NN = 25;
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .8;

bvec = 6:20;
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

% Plot of error with MFS
subplot(1,2,1)
loglog(4*(bvec-1),errMFS)
ylabel('error')
xlabel('Collocation points')
title('Solution via MFS')
set(gca,'xtick',4*(bvec-1))

% Solve the system with GaussQRr with partial regression
fsol = @(x,y) exp(.5*x+.5).*cos(pi/4*(y+1));
f = @(x,y) (1/4-pi^2/16)*fsol(x,y);

ptsEVAL = pick2Dpoints([-1 -1],[1 1],NN);
usol = fsol(ptsEVAL(:,1),ptsEVAL(:,2));

alpha = 1;
ep = 1e-9;

m = 1;
errvecR2D = [];
Nvec = [];
for N=5:10
    % Choose some possible collocation points and eliminate duplicates
    x = [pick2Dpoints([-1 -1],[1 1],N,'cheb');pick2Dpoints([-1 -1],[1 1],N,'halton')];
    x = unique(1e-8*ceil(1e8*x),'rows');
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
    GQR.N = length(x);
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

pause


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here now is the second example
% This time we use the same function as before, but the domain
% has the top quadrant missing, making it a nonconvex domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fsol = @(x,y) exp(x).*cos(y);
Lfs = @(r) -1/(2*pi)*log(r);

bvec = 14:6:50;
bNvec = [];
errMFS = [];

m = 1;
for bN=bvec
    % Choose the collocation points
    x = [pick2Dpoints([-1 -1],[1 1],bN);pick2Dpoints([0 0],[1 1],ceil(bN/2))];
    bx = find( x(:,1)==-1 | x(:,2)==-1 | (x(:,1)==1 & x(:,2)<=0) | (x(:,2)==1 & x(:,1)<=0) | (x(:,1)>=0 & x(:,2)==0) | (x(:,2)>=0 & x(:,1)==0) );
    ptsMFScoll = unique(1e-8*ceil(1e8*x(bx,:)),'rows')*diag([.5 pi/4]) + ones(size(bx,1),1)*[.5 pi/4];
    bNvec(m) = size(ptsMFScoll,1);
    
    ptsMFSsource = 3*[cos(linspace(-pi,pi,bNvec(m)))',sin(linspace(-pi,pi,bNvec(m)))'] + ones(bNvec(m),1)*[.5,pi/4];
    % Find some sample points to evaluate the error at
    x = pick2Dpoints([-1 -1],[1 1],NN);
    bx = find( x(:,1)<=0 | x(:,2)<=0 );
    ptsEVAL = x(bx,:)*diag([.5 pi/4]) + ones(size(bx,1),1)*[.5 pi/4];
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

% Plot of error with MFS
subplot(1,2,1)
loglog(bNvec,errMFS)
ylabel('error')
xlabel('Collocation points')
title('Solution via MFS')
set(gca,'xtick',bNvec)
xlim([min(bNvec),max(bNvec)])


% GaussQR solution
fsol = @(x,y) exp(.5*x+.5).*cos(pi/4*(y+1));
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
    % Determine the collocation points for the L-shaped domain
    possible_x = [pick2Dpoints([-1 -1],[1 1],N,'cheb');pick2Dpoints([-1 -1],[1 1],N,'halton');pick2Dpoints([0 0],[1 1],ceil(N/2))];
    % Get rid of duplicate points, which can happen
    x = unique(1e-8*ceil(1e8*possible_x),'rows');
    bx = find( x(:,1)<=0 | x(:,2)<=0 );
    x = x(bx,:);
    % Determine which points are boundary points, and which are interior
    b = find( abs(x(:,1))==1 | abs(x(:,2))==1 | (x(:,1)==0 & x(:,2)>=0) | (x(:,1)>=0 & x(:,2)==0));
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
    GQR.N = Nvec(m);
    GQR.coef = coef;

    uGQR = rbfqr_eval(GQR,ptsEVAL);
    errvecR2D(m) = errcompute(uGQR,usol);
    m = m + 1;
end

subplot(1,2,2)
loglog(Nvec,errvecR2D)
ylabel('error')
xlabel('Collocation points')
title('Solution via GaussQR')
xlim([min(Nvec),max(Nvec)])
set(gca,'xtick',Nvec)




pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here now is the third example
% This time we use the same function and domain as before,
% but the BC are mixed Dirichlet and Neumann
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fsol = @(x,y) exp(x).*cos(y);
fdysol = @(x,y) -exp(x).*sin(y);
Lfs = @(r) -1/(2*pi)*log(r);
Ldyfs = @(x,y) -1./(2*pi*DistanceMatrix(x,y).^2).*abs(DifferenceMatrix(x(:,2),y(:,2)));

bvec = 20:10:90;
bNvec = [];
errMFS = [];

m = 1;
for bN=bvec
    % Choose the collocation points
    x = [pick2Dpoints([-1 -1],[1 1],bN);pick2Dpoints([0 0],[1 1],ceil(bN/2))];
    bx = find( x(:,1)==-1 | x(:,2)==-1 | (x(:,1)==1 & x(:,2)<=0) | (x(:,2)==1 & x(:,1)<=0) | (x(:,1)>=0 & x(:,2)==0) | (x(:,2)>=0 & x(:,1)==0) );
    ptsMFScoll = unique(1e-8*ceil(1e8*x(bx,:)),'rows')*diag([.5 pi/4]) + ones(size(bx,1),1)*[.5 pi/4];
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

% Plot of error with MFS
subplot(1,2,1)
loglog(bNvec,errMFS)
ylabel('error')
xlabel('Collocation points')
title('Solution via MFS')
set(gca,'xtick',bNvec)
xlim([min(bNvec),max(bNvec)])


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
    [ep,alpha,Marr] = rbfsolveprep(1,x,ep,alpha);
    % Build the collocation matrix
    phiMat2d = rbfphi(Marr,x(bi,:),ep,alpha,[2,0])+rbfphi(Marr,x(bi,:),ep,alpha,[0,2]);
    phiMat = rbfphi(Marr,x(bD,:),ep,alpha);
    phiMat1dy = rbfphi(Marr,x(bN,:),ep,alpha,[0,1]);
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

    uGQR = rbfqr_eval(GQR,ptsEVAL);
    errvecR2D(m) = errcompute(uGQR,usol);
    m = m + 1;
end

subplot(1,2,2)
loglog(Nvec,errvecR2D)
ylabel('error')
xlabel('Collocation points')
title('Solution via GaussQR')
xlim([min(Nvec),max(Nvec)])
set(gca,'xtick',Nvec)

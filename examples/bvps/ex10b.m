% ex10b.m
% Incorporating the stable HS-SVD basis into MPS
%    This problem takes place on an L-shaped domain
% In honor of Graeme Fairweather's 70th birthday
% This appears in:
%    Using Gaussian eigenfunctions to solve boundary value problems
%    M. McCourt, Advances in Applied Mathematics and Mechanics, 5:569-594, 2013.
%
% The problem we look at is 
%     Lap(u) - lambda^2*u = f   -on- interior
%     u = g                     -on- boundary
%   solution: u(x,y) = sin(x^2+y)
%   domain: L shaped region (-1<x<0 & -1<y<1)+(-1<x<1 & -1<y<0)
%   fictitious boundary: 1.1 times real boundary

global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

% This is the wavenumber (maybe) for the Helmholtz problem
lambda = 3;

% This is the true solution of the problem
fsol = @(x,y) sin(x.^2+y);
fsoldy = @(x,y) cos(x.^2+y);
% The source associated with the true solution
f = @(x,y) 2*cos(x.^2+y)-4*x.^2.*sin(x.^2+y)-sin(x.^2+y)-lambda^2*fsol(x,y);

% fsol = @(x,y) exp(x+y);
% f = @(x,y) (2-lambda^2)*exp(x+y);
% The fundamental solution of the Helmholtz problem
Hfs = @(r) besselk(0,lambda*r)/(2*pi);
Hdyfs = @(x,y) -lambda*DifferenceMatrix(x(:,2),y(:,2)).*besselk(1,lambda*DistanceMatrix(x,y))./DistanceMatrix(x,y);

% The number of error evalution points in each dimension
NN = 30;

% The length of the GaussQRr regression
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .5;
% The global scale parameter for GaussQR
alpha = 1;
% The shape parameter commonly associated with RBFs
ep = 1e-6;

bvec = 10:5:70;
errMPS = zeros(size(bvec));
errMPSimp = zeros(size(bvec));
errGQR = zeros(size(bvec));
bNvec = zeros(size(bvec));
Nvec = zeros(size(bvec));

warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:rankDeficientMatrix')

h_waitbar = waitbar(0,'Initializing');

m = 1;
for bnum=bvec
    % Choose the collocation points
    x = [pick2Dpoints([-1 -1],[1 1],bnum);pick2Dpoints([0 0],[1 1],ceil(bnum/2))];
    bx = ( x(:,1)==-1 | x(:,2)==-1 | (x(:,1)==1 & x(:,2)<=0) | (x(:,2)==1 & x(:,1)<=0) | (x(:,1)>=0 & x(:,2)==0) | (x(:,2)>=0 & x(:,1)==0) );
    x = unique(1e-8*ceil(1e8*x(bx,:)),'rows');
    bNvec(m) = size(x,1);
    bN = find( x(:,2)==-1 );
    bD = setdiff(1:bNvec(m),bN);
    ptsMFScoll_D = x(bD,:);
    ptsMFScoll_N = x(bN,:);
    ptsMFScoll = [ptsMFScoll_D;ptsMFScoll_N];
    N = ceil(sqrt(bNvec(m)));
    
    % Find the MFS source points
    ptsMFSsource = 1/(lambda^2+N)*(ptsMFScoll.*(abs(ptsMFScoll)==1) + (ptsMFScoll==0)) + ptsMFScoll;
%     tv = linspace(0,2*pi,9)';
%     rv = [l2,1+l2,sqrt(2)+l2,1+l2,sqrt(2)+l2,1+l2,sqrt(2)+l2,1+l2,l2]';
%     tt = linspace(0,2*pi,bNvec(m))';
%     rr = interp1(tv,rv,tt);
%     ptsMFSsource = [rr.*cos(tt+pi/4),rr.*sin(tt+pi/4)];
    
    % Find some sample points to evaluate the error at
    x = pick2Dpoints([-1 -1],[1 1],NN);
    bx = ( x(:,1)<=0 | x(:,2)<=0 );
    ptsEVAL = x(bx,:);
    % Evaluate the true solution at the sample points
    usol = fsol(ptsEVAL(:,1),ptsEVAL(:,2));
    
    % Choose the GaussQR interpolation points
    possible_x = [pick2Dpoints([-1 -1],[1 1],N,'cheb');pick2Dpoints([-1 -1],[1 1],N,'halton')];
    % Get rid of duplicate points, which can happen
    x = unique(1e-8*ceil(1e8*possible_x),'rows');
    bx = find( x(:,1)<=0 | x(:,2)<=0 );
    x = x(bx,:);
    % Determine which points are boundary points, and which are interior
    b = find( abs(x(:,1))==1 | abs(x(:,2))==1 | (x(:,1)==0 & x(:,2)>=0) | (x(:,1)>=0 & x(:,2)==0) );
    bi = setdiff(1:size(x,1),b)';
    ptsGQR = x(bi,:);
    Nvec(m) = size(ptsGQR,1);
    
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
    
    % Store the approximation length to be used later
    M_MPS = size(GQR.Marr,2);
    
    % Now enforce the boundary with MFS
    uPonBDY_D = gqr_eval(GQR,ptsMFScoll_D);
    uPonBDY_N = gqr_eval(GQR,ptsMFScoll_N,[0,1]);
    rhs_D = fsol(ptsMFScoll_D(:,1),ptsMFScoll_D(:,2)) - uPonBDY_D;
    rhs_N = fsoldy(ptsMFScoll_N(:,1),ptsMFScoll_N(:,2)) - uPonBDY_N;
    A_coll_D = Hfs(DistanceMatrix(ptsMFScoll_D,ptsMFSsource));
    A_coll_N = Hdyfs(ptsMFScoll_N,ptsMFSsource);
    coefMFS = [A_coll_D;A_coll_N]\[rhs_D;rhs_N];
    
    % Consider the problem with just GaussQR for comparison
    % Different boundary points need to be used here than MFS
    N = floor(.3*N^2/4); % Number of points on each side
    x = [ [pickpoints(-1,1,N,'cheb'),-ones(N,1)];[-ones(N,1),pickpoints(-1,1,N,'cheb')];[pickpoints(-1,1,N,'cheb'),ones(N,1)];[ones(N,1),pickpoints(-1,1,N,'cheb')] ];
%     x = pick2Dpoints([-1,-1],[1 1],6);
    x = unique(1e-8*ceil(1e8*x),'rows');
    x = x(any(abs(x)==1,2),:);
    ib = find((x(:,1)>0) & (x(:,2)==1));
    x(ib,2) = zeros(size(x(ib,2)));
    ib = find((x(:,1)==1) & (x(:,2)>0));
    x(ib,1) = zeros(size(x(ib,1)));
    bN = find( x(:,2)==-1 );
    bD = setdiff(1:size(x,1),bN);
    ptsBDY_D = x(bD,:);
    ptsBDY_N = x(bN,:);
    ptsBDY = x;
    
    % Solve the system using these different boundary points
    ptsFULL = [ptsGQR;ptsBDY];
    [ep,alpha,Marr] = gqr_solveprep(1,ptsFULL,ep,alpha,M_MPS);
    phiMat = gqr_phi(Marr,ptsGQR,ep,alpha);
    phiMatBC_D = gqr_phi(Marr,ptsBDY_D,ep,alpha);
    phiMatBC_N = gqr_phi(Marr,ptsBDY_N,ep,alpha,[0,1]);
    phiMat2d = gqr_phi(Marr,ptsGQR,ep,alpha,[2,0])+gqr_phi(Marr,ptsGQR,ep,alpha,[0,2]);
    A = [phiMat2d - lambda^2*phiMat;phiMatBC_N;phiMatBC_D];
    rhs_L = f(ptsGQR(:,1),ptsGQR(:,2));
    rhs_N = fsoldy(ptsBDY_N(:,1),ptsBDY_N(:,2));
    rhs_D = fsol(ptsBDY_D(:,1),ptsBDY_D(:,2));
    rhs = [rhs_L;rhs_N;rhs_D];
    coef = A\rhs;
    
    % Fill the GQR object with all the values it needs
    GQRfull.reg = 1;
    GQRfull.Marr = Marr;
    GQRfull.alpha = alpha;
    GQRfull.ep = ep;
    GQRfull.N = size(ptsFULL,1);
    GQRfull.coef = coef;
    
    % Also consider improved MFS with full GaussQR
    uPonBDY_D = gqr_eval(GQRfull,ptsMFScoll_D);
    uPonBDY_N = gqr_eval(GQRfull,ptsMFScoll_N,[0,1]);
    rhs_D = fsol(ptsMFScoll_D(:,1),ptsMFScoll_D(:,2)) - uPonBDY_D;
    rhs_N = fsoldy(ptsMFScoll_N(:,1),ptsMFScoll_N(:,2)) - uPonBDY_N;
    A_coll_D = Hfs(DistanceMatrix(ptsMFScoll_D,ptsMFSsource));
    A_coll_N = Hdyfs(ptsMFScoll_N,ptsMFSsource);
    coefMFSimp = [A_coll_D;A_coll_N]\[rhs_D;rhs_N];
    
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
    
    waitbar(m/length(bvec),h_waitbar,sprintf('num bdy points: %d',bnum))
    m = m+1;
end

waitbar(1,h_waitbar,'Plotting')

warning('on','MATLAB:nearlySingularMatrix')
warning('on','MATLAB:rankDeficientMatrix')

% Plot of error with MFS
h = figure;
loglog(bNvec,[errGQR;errMPS;errMPSimp],'linewidth',3)
ylabel('Absolute 2-norm error')
xlabel('Number of interior points')
xlim([min(bNvec),max(bNvec)])
ylim([1e-10,1e2])
legend('HS-SVD','MPS','MPS+HS-SVD','location','southwest')
set(gca,'xtick',[36,50,80,140,250])
set(gca,'ytick',[1e-10,1e-7,1e-4,1e-1,1e2])

close(h_waitbar)
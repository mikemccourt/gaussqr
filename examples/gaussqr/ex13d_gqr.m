% ex13d_gqr.m
% Incorporating GaussQR into MPS
% In honor of Graeme Fairweather's 70th birthday
% The first problem we'll be looking at is 
%     Lap(u) - lambda^2*u = f   -on- interior
%     u = g                     -on- boundary
%   solution: u(x,y) = sin(x^2+y)
%   domain: L shaped region (-1<x<0 & -1<y<1)+(-1<x<1 & -1<y<0)
%   fictitious boundary: 1.1 times real boundary

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

% This is the wavenumber (maybe) for the Helmholtz problem
lambda = 3;

% This is the true solution of the problem
fsol = @(x,y) sin(x.^2+y);
% The source associated with the true solution
f = @(x,y) 2*cos(x.^2+y)-4*x.^2.*sin(x.^2+y)-sin(x.^2+y)-lambda^2*fsol(x,y);

% fsol = @(x,y) exp(x+y);
% f = @(x,y) (2-lambda^2)*exp(x+y);
% The fundamental solution of the Helmholtz problem
Hfs = @(r) besselk(0,lambda*r)/(2*pi);

% The number of error evalution points in each dimension
NN = 35;

% The length of the GaussQRr regression
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .5;
% The global scale parameter for GaussQR
alpha = 1;
% The shape parameter commonly associated with RBFs
ep = 1e-5;

bvec = 10:5:70;
errMPS = [];
errMPSimp = [];
errGQR = [];
bNvec = [];

m = 1;
for bN=bvec
    % Choose the collocation points
    x = [pick2Dpoints([-1 -1],[1 1],bN);pick2Dpoints([0 0],[1 1],ceil(bN/2))];
    bx = find( x(:,1)==-1 | x(:,2)==-1 | (x(:,1)==1 & x(:,2)<=0) | (x(:,2)==1 & x(:,1)<=0) | (x(:,1)>=0 & x(:,2)==0) | (x(:,2)>=0 & x(:,1)==0) );
    ptsMFScoll = unique(1e-8*ceil(1e8*x(bx,:)),'rows');
    bNvec(m) = size(ptsMFScoll,1);
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
    bx = find( x(:,1)<=0 | x(:,2)<=0 );
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
    b = find( abs(x(:,1))==1 | abs(x(:,2))==1 | (x(:,1)==0 & x(:,2)>=0) | (x(:,1)>=0 & x(:,2)==0));
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
    
    % Now enforce the boundary with MFS
    uPonBDY = gqr_eval(GQR,ptsMFScoll);
    rhs = fsol(ptsMFScoll(:,1),ptsMFScoll(:,2)) - uPonBDY;
    DM_coll = DistanceMatrix(ptsMFScoll,ptsMFSsource);
    A_coll = Hfs(DM_coll);
    coefMFS = A_coll\rhs;
    
    % Evaluate the error
    uP_eval = gqr_eval(GQR,ptsEVAL);
    DM_eval = DistanceMatrix(ptsEVAL,ptsMFSsource);
    A_eval = Hfs(DM_eval);
    uF_eval = A_eval*coefMFS;
    errMPS(m) = errcompute(uF_eval+uP_eval,usol);
    
    % Consider the problem with just GaussQR for comparison
    % Different boundary points need to be used here than MFS
    N = floor(.4*N^2/4); % Number of points on each side
    x = [ [pickpoints(-1,1,N,'cheb'),-ones(N,1)];[-ones(N,1),pickpoints(-1,1,N,'cheb')];[pickpoints(-1,1,N,'cheb'),ones(N,1)];[ones(N,1),pickpoints(-1,1,N,'cheb')] ];
    x = unique(1e-8*ceil(1e8*x),'rows');
    ib = find((x(:,1)>0) & (x(:,2)==1));
    x(ib,2) = zeros(size(x(ib,2)));
    ib = find((x(:,1)==1) & (x(:,2)>0));
    x(ib,1) = zeros(size(x(ib,1)));
    ptsBDY = x;
    
    % Solve the system using these different boundary points
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
    
    fprintf('%d ',m)
    m = m+1;
end

% Plot of error with MFS
loglog(bNvec,[errMPS;errGQR;errMPSimp],'linewidth',3)
ylabel('RMS error')
xlabel('Collocation points')
xlim([min(bNvec),max(bNvec)])
legend('MPS','GaussQR','MPS+GaussQR','location','southwest')
%set(gca,'xtick',bNvec)

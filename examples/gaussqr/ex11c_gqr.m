% ex11
% This example is meshfree coupling between domains
% The content is similar to ex11, but here I'm trying to
% make a convergence picture
%
% The problem we are considering is
%   u_t - u_xx = f
% with Dirichlet BC and second order interface conditions
% f is determined by the solution given
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 4;
regfrac = GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC;

AN = [4,8,16,32,64,128];
BN = [4,8,16,32,64,128];
epopt1 = [.989,.9886,1.0512,.9267,.8861,.8030];
rfrac1 = [.9,.5,.3,.3,.2,.1];
ind = 6; % Handles the number of points
count = 1; % Handles the number of time steps
restart = 30; % GMRES restart value

% Number of points included in the coupling region
couplewidth = 1;
% Technique used for coupling
% 1 - Traditional finite difference
% 2 - Meshfree
% 3 - Both
Fcstyle = 3;

% This is needed for RBF direct
% I guess it isn't really needed, but I'll keep it for now
% This is the Gaussian and its x derivative
rbf = @(ep,r) exp(-(ep*r).^2);
rbfdx = @(ep,r,dx) -2*ep^2*dx.*exp(-(ep*r).^2);
% The Gaussian shape parameter
epv = 1e-3;
% The global scale parameter for GaussQR
alpha = 1;

% This allows you to test multiple epsilon values
%  and see which is the best for error
% The number of epsilon points to consider
FM = 21;
epvec = [];
% Uncomment this line to use fminbnd to find an 
%  optimal choice of epsilon using those values
%  as the bounds of a search interval
%epvec = [.5,1.3];
% Uncomment this line as well to evaluate the error
%  all the epsilons in epvec.  The results are stored
%  in the vector Fapx
%epvec = logspace(-.5,.5,FM);

% True solution and source
testopt = 1;
switch testopt
    case 1
        testname = 'polynomial';
        yf = @(x,t) exp(-t)*(1-0.5*(x(:,1).^2+x(:,2).^2));
        ff = @(x,t) exp(-t)*(1+0.5*(x(:,1).^2+x(:,2).^2));
    case 2
        testname = 'sinh';
        yf = @(x,t) exp(-t)*sinh(x(:,1)+x(:,2));
        ff = @(x,t) -3*yf(x,t);
end

% I may turn this back on eventually
warning off

% Setting up the boundary value problem discretization
dt = 0.01/(2^(count-1));
nsteps = 2*(2^(count-1));
FDerrt = zeros(1,nsteps);
MFerrt = zeros(1,nsteps);

% This is the 4th order finite difference stencil
stencil9 = [1/6 2/3 1/6  2/3 -10/3 2/3  1/6 2/3 1/6];
% A small fudge factor to allow for fuzzy math
%  in computing the grid
couplebuffer = 0.1/(2^(ind-1));

% Set up the grid of the problem
Ald = [-1 0];Aud = [0 1];
Bld = [0 0];Bud = [1 1];
AM = AN(ind);
BM = BN(ind);

% Just an evenly spaced grid, for now
Ax = pick2Dpoints(Ald,Aud,AN(ind)*ones(2,1));
Bx = pick2Dpoints(Bld,Bud,BN(ind)*ones(2,1));
AMM = size(Ax,1);
BMM = size(Bx,1);
Adeltax = 1/(AM-1);
Bdeltax = 1/(BM-1);

% Need to find the points which only A needs, the points in the coupling
% region, and the points on the interface between models
% We also store the interior points, which are the points governed by the
% PDE so that we can map the stencil over
Aifapts = find(abs(Ax(:,1)-0)<couplebuffer)';
Acoupts = setdiff(find(abs(Ax(:,1)-0)<=Adeltax*(couplewidth+couplebuffer)),Aifapts);
Aonlpts = setdiff(1:AMM,[Aifapts,Acoupts]);
Aintpts = find(1-any((Ax==1)+(Ax==0)+(Ax==-1),2))';

% Need to set up the data in an appropriate multiphysics way
% We'll use Aonly,Acoupling,Ainterface
% Given the reordering, we'll need to make an adjustment in the finite
% difference stencil as well
Ashuffle = [Aonlpts,Acoupts,Aifapts];
[temp,Aunshuffle] = sort(Ashuffle);
Ax = Ax(Ashuffle,:);
Astenoff9 = [-AM-1 -AM -AM+1  -1 0 1  AM-1 AM AM+1];
Asten = zeros(length(Aintpts),length(Astenoff9));
for k=1:length(Aintpts)
  Asten(k,:) = Aunshuffle(Astenoff9+Aintpts(k));
end

% Refind the sets of points, and make sure they are sorted
% The sorting will automatically occur because of the find function
Aifapts = find(abs(Ax(:,1)-0)<couplebuffer)';
Acoupts = setdiff(find(abs(Ax(:,1)-0)<=Adeltax*(couplewidth+couplebuffer)),Aifapts);
Aonlpts = setdiff(1:AMM,[Aifapts,Acoupts]);
Aintpts = find(1-any((Ax==1)+(Ax==0)+(Ax==-1),2))';
Aboupts = setdiff(1:AMM,Aintpts);

% Have to account for the reordering because the find function used to
% generate Aintpts will return a sorted list.  Asten must be sorted too.
[temp,order] = sort(Asten(:,ceil(length(Astenoff9)/2)));
Asten = Asten(order,:);


% Do all that same stuff for the B model
Bifapts = find(abs(Bx(:,1)-0)<couplebuffer)';
Bcoupts = setdiff(find(abs(Bx(:,1)-0)<=Bdeltax*(couplewidth+couplebuffer)),Bifapts);
Bonlpts = setdiff(1:BMM,[Bifapts,Bcoupts]);
Bintpts = find(1-any((Bx==1)+(Bx==0)+(Bx==-1),2))';

Bshuffle = [Bonlpts,Bcoupts,Bifapts];
[temp,Bunshuffle] = sort(Bshuffle);
Bx = Bx(Bshuffle,:);
Bstenoff9 = [-BM-1 -BM -BM+1  -1 0 1  BM-1 BM BM+1];
Bsten = zeros(length(Bintpts),length(Bstenoff9));
for k=1:length(Bintpts)
  Bsten(k,:) = Bunshuffle(Bstenoff9+Bintpts(k));
end

Bifapts = find(abs(Bx(:,1)-0)<couplebuffer)';
Bcoupts = setdiff(find(abs(Bx(:,1)-0)<=Bdeltax*(couplewidth+couplebuffer)),Bifapts);
Bonlpts = setdiff(1:BMM,[Bifapts,Bcoupts]);
Bintpts = find(1-any((Bx==1)+(Bx==0)+(Bx==-1),2))';
Bboupts = setdiff(1:BMM,Bintpts);

[temp,order] = sort(Bsten(:,ceil(length(Bstenoff9)/2)));
Bsten = Bsten(order,:);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we need to form the coupled system, under the assumption
%  that both models can pass info back and forth as needed
%
% We need to determine what the ordering is going to be
%  [Aonly,Bonly,Acoupling,Ainterface,Bcoupling,Binterface]
% The shuffle described below would occur if the matrix were ordered
%  Fmat = [Amat,zeros(BMM,AMM);zeros(AMM,BMM),Bmat];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fshuffle = [Aonlpts,Bonlpts+AMM,Acoupts,Aifapts,Bcoupts+AMM,Bifapts+AMM];
[temp,Funshuffle] = sort(Fshuffle);
k = 0;
FAonl = k+1:k+length(Aonlpts);
k = k + length(Aonlpts);
FBonl = k+1:k+length(Bonlpts);
k = k + length(Bonlpts);
FAcou = k+1:k+length(Acoupts);
k = k + length(Acoupts);
FAifa = k+1:k+length(Aifapts);
k = k + length(Aifapts);
FBcou = k+1:k+length(Bcoupts);
k = k + length(Bcoupts);
FBifa = k+1:k+length(Bifapts);

% Set up the points needed for the problem
Fx = [Ax;Bx];
Fx = Fx(Fshuffle,:);
FMM = AMM + BMM;

% Set the initial conditions, provided by the true solution
Aunew = yf(Ax,0);
Bunew = yf(Bx,0);

for tn=1:nsteps
    t = tn*dt;
    Fsol  = yf(Fx,t); % true solution, for computing error

    Amat = 1.0/dt*speye(AMM);
    Arhs = zeros(AMM,1);
    Arf = ff(Ax,t);
    Aubc = yf(Ax,t);
    Auold = Aunew;
    m = 1;
    for k=Aintpts
        Amat(k,Asten(m,:)) = Amat(k,Asten(m,:)) - stencil9/(Adeltax^2);
        Arhs(k) = Arf(k) + 1/dt*Auold(k);
        m = m + 1;
    end
    for k=Aboupts
        Amat(k,k) = 1.0;
        Arhs(k) = Aubc(k);
    end

    Bmat = 1/dt*speye(BMM);
    Brhs = zeros(BMM,1);
    Brf = ff(Bx,t);
    Bubc = yf(Bx,t);
    Buold = Bunew;
    m = 1;
    for k=Bintpts
        Bmat(k,Bsten(m,:)) = Bmat(k,Bsten(m,:)) - stencil9/(Bdeltax^2);
        Brhs(k) = Brf(k) + 1/dt*Buold(k);
        m = m + 1;
    end
    for k=Bboupts
        Bmat(k,k) = 1.0;
        Brhs(k) = Bubc(k);
    end

    % Now that we have the pieces, we can form the fully coupled system
    Fmat = [Amat,zeros(BMM,AMM);zeros(AMM,BMM),Bmat];
    Fmat = Fmat(Fshuffle,Fshuffle);
    Frhs = [Arhs;Brhs];
    Frhs = Frhs(Fshuffle);
    
    % Make the replacements associated with the coupling
    % Then solve the system and compute the error
    if Fcstyle==1 | Fcstyle==3
        % Copy the existing non-coupled matrix
        FmatFD = Fmat;

        FmatFD(FAifa,:) = zeros(size(FmatFD(FAifa,:)));
        for k=1:length(FAifa)
            FmatFD(FAifa(k),[FAifa(k),FBifa(k)]) = [1 -1];
        end
        Frhs(FAifa) = zeros(size(FAifa));
        
        FmatFD(FBifa,:) = zeros(size(FmatFD(FBifa,:)));
        for k=1:length(FBifa)
            if couplewidth==1
                FmatFD(FBifa(k),[FAcou(k),FAifa(k),FBcou(k),FBifa(k)]) = [[-1 1]/Adeltax [-1 1]/Bdeltax];
            elseif couplewidth==2
                FmatFD(FBifa(k),[FAcou(2*k-[1,0]),FAifa(k),FBcou(2*k-[0,1]),FBifa(k)]) = [[-1 4 -3]/(2*Adeltax) [-1 4 -3]/(2*Bdeltax)];
            else
                error('couplewidth should be 1 or 2')
            end
        end
        Frhs(FBifa) = zeros(size(FBifa));

        % Solve the system and compute the error on the grid
        Funew = FmatFD\Frhs;
        FDerrt(tn) = errcompute(Funew,Fsol);
    end
    if Fcstyle==2 | Fcstyle==3 % Meshfree coupling
        % Copy the existing non-coupled matrix
        FmatMF = Fmat;
        if not(exist('FmatFD'))
            FmatFD = [];
        end

        % This couples the values at the interface
        % Need to work on mismatched mesh eventually
        FmatMF(FAifa,:) = zeros(size(FmatMF(FAifa,:)));
        for k=1:length(FAifa)
            FmatMF(FAifa(k),[FAifa(k),FBifa(k)]) = [1 -1];
        end
        Frhs(FAifa) = zeros(size(FAifa));

        AMFpts = Fx([FAcou,FAifa],:);
        BMFpts = Fx([FBcou,FBifa],:);
        ap = length(AMFpts);
        
        % This section couples the derivatives at the interface
        % If requested, this will do a parameter search
        % Otherwise, skip this for loop
        if not(isempty(epvec))
            Marr = rbfformMarr([couplewidth+1;ap],[],floor(regfrac*ap));
            epstruct.Marr = Marr;
            epstruct.alpha = alpha;
            epstruct.AMFpts = AMFpts;
            epstruct.BMFpts = BMFpts;
            epstruct.FmatMF = FmatMF;
            epstruct.FAcou = FAcou;
            epstruct.FAifa = FAifa;
            epstruct.FBcou = FBcou;
            epstruct.FBifa = FBifa;
            epstruct.Frhs = Frhs;
            epstruct.restart = restart;
            epstruct.Fsol = Fsol;
            epstruct.Fx = Fx;
            epstruct.FmatFD = FmatFD;

            epv = fminbnd(@(ep)ex11cTestFunc(ep,epstruct),min(epvec),max(epvec));
        end

        % This will loop over all the values in epvec and
        % compute the error to let you know which is the best
        Fapx = zeros(FM,1);
        k = 1;
        if length(epvec)>2
            for ep=epvec
                Marr = rbfformMarr([couplewidth+1;ap],[],floor(regfrac*ap));
                phiA = rbfphi(Marr,AMFpts,ep,alpha);
                phiAx = rbfphi(Marr,Fx(FBifa,:),ep,alpha,[1 0]);
                phiB = rbfphi(Marr,BMFpts,ep,alpha);
                phiBx = rbfphi(Marr,Fx(FBifa,:),ep,alpha,[1 0]);
    
                FmatMF(FBifa,:) = zeros(size(FmatMF(FBifa,:)));
                FmatMF(FBifa,[FAcou,FAifa]) = phiAx/phiA;
                FmatMF(FBifa,[FBcou,FBifa]) = -phiBx/phiB;
                Frhs(FBifa) = zeros(size(FBifa));
                
                [Funew,cnvg] = gmres(FmatMF,Frhs,restart,[],[],FmatFD);
                Fapx(k) = errcompute(Funew,Fsol);

%                fprintf('%d\t%g\t%g\n',k,ep,Fapx(k))
                k = k+1;
            end
            [err,epi] = min(Fapx);
            epv = epvec(epi);
        end

        Marr = rbfformMarr([couplewidth+1;ap],[],floor(regfrac*ap));
        phiA = rbfphi(Marr,AMFpts,epv,alpha);
        phiAx = rbfphi(Marr,Fx(FBifa,:),epv,alpha,[1 0]);
        phiB = rbfphi(Marr,BMFpts,epv,alpha);
        phiBx = rbfphi(Marr,Fx(FBifa,:),epv,alpha,[1 0]);

        FmatMF(FBifa,:) = zeros(size(FmatMF(FBifa,:)));
        FmatMF(FBifa,[FAcou,FAifa]) = phiAx/phiA;
        FmatMF(FBifa,[FBcou,FBifa]) = -phiBx/phiB;
        Frhs(FBifa) = zeros(size(FBifa));

        % Solve the system and compute the error on the grid
        [Funew,cnvg] = gmres(FmatMF,Frhs,restart,[],[],FmatFD);
        MFerrt(tn) = errcompute(Funew,Fsol);
    end
   
    % Unscatter the components back to A and B
    temp = Funew(Funshuffle);
    Aunew = temp(1:AMM);
    Bunew = temp(AMM+1:end);
end
warning on




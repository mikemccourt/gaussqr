% ex11
% This example is meshfree coupling between domains
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2; % Absolute error

epvecd = logspace(-1,1,20);
epvecr = logspace(-1,1,20);
AN = [4,8,16,32,64,128];
BN = [4,8,16,32,64,128];
rbf = @(ep,x) exp(-(ep*x).^2);
ep = 1.0;

if not(exist('ind')) % Consider the first index (looping only)
  ind = 1;
end

if not(exist('count')) % Work with the basic time step
  count = 1;
end

dt = 0.01/(2^(count-1));
nsteps = 2*(2^(count-1));
stencil9 = [1/6 2/3 1/6  2/3 -10/3 2/3  1/6 2/3 1/6];
errt = zeros(1,nsteps);
couplebuffer = 0.1/(2^(ind-1)); % accounts for fuzzy math
couplewidth = 2; % How many points included in coupling

yf = @(x,t) exp(-t)*(1-0.5*(x(:,1).^2+x(:,2).^2));
ff = @(x,t) exp(-t)*(1+0.5*(x(:,1).^2+x(:,2).^2));

Ald = [-1 0];Aud = [0 1];
Bld = [0 0];Bud = [1 1];
AM = AN(ind);
BM = BN(ind);

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

% Set the initial conditions, provided by the true solution
Aunew = yf(Ax,0);
Bunew = yf(Bx,0);
Aerrt = zeros(1,nsteps);
Berrt = zeros(1,nsteps);

% Step through time
for tn=1:nsteps
  t = tn*dt;

  % Initialize the values of the matrix and RHS
  Amat = 1/dt*speye(AMM);
  Arhs = zeros(AMM,1);
  
  % This is the source term evaluated on the grid
  Arf = ff(Ax,t);
  % The boundary conditions are derived from the true solution
  Aubc = yf(Ax,t);
  % Store the old time step for use below
  Auold = Aunew;

  % First consider only the PDE part
  m = 1;
  for k=Aintpts
      Amat(k,Asten(m,:)) = Amat(k,Asten(m,:)) - stencil9/(Adeltax^2);
      Arhs(k) = Arf(k) + 1/dt*Auold(k);
      m = m + 1;
  end
  % Now consider the boundary condition part
  for k=Aboupts
      Amat(k,k) = 1.0;
      Arhs(k) = Aubc(k);
  end
  
  % Solve the system and compute the error on the grid
  Aunew = Amat\Arhs;
  Aerrt(tn) = errcompute(Aunew,Aubc);


  % Do the same thing with the B system
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

  Bunew = Bmat\Brhs;
  Berrt(tn) = errcompute(Bunew,Bubc);
end

% Now that we've considered the result of solving the systems individually,
% we want to determine what effect is caused by solving them with a
% specific coupling strategy.

% We need to determine what the ordering is going to be
% [Aonly,Bonly,Acoupling,Ainterface,Bcoupling,Binterface]
% The shuffle described below would occur if the matrix were ordered
% Fmat = [Amat,zeros(BMM,AMM);zeros(AMM,BMM),Bmat];
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

% The simplest coupling style just matches values and derivatives
Fcstyle = 1;
Ferrt = zeros(1,nsteps);

Fx = [Ax;Bx];
Fx = Fx(Fshuffle,:);
Aunew = yf(Ax,0);
Bunew = yf(Bx,0);
FMM = AMM + BMM;

for tn=1:nsteps
    t = tn*dt;
    Fsol  = yf(Fx,t);

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
    if Fcstyle==1
        Fmat(FAifa,:) = zeros(size(Fmat(FAifa,:)));
        for k=1:length(FAifa)
            Fmat(FAifa(k),[FAifa(k),FBifa(k)]) = [1 -1];
            Frhs(FAifa(k)) = 0;
        end
        Fmat(FBifa,:) = zeros(size(Fmat(FBifa,:)));
        for k=1:length(FBifa)
            Fmat(FBifa(k),[FAcou(k),FAifa(k),FBcou(k),FBifa(k)]) = [[-1 1]/Adeltax [-1 1]/Bdeltax];
            Frhs(FBifa(k)) = 0;
        end
    end
    
    % Solve the system and compute the error on the grid
    Funew = Fmat\Frhs;
    Ferrt(tn) = errcompute(Funew,Fsol);
    
    % Unscatter the components back to A and B
    temp = Funew(Funshuffle);
    Aunew = temp(1:AMM);
    Bunew = temp(AMM+1:end);
end

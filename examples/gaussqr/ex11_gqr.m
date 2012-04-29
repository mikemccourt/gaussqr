% ex11
% This example is meshfree coupling between domains
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2; % Absolute error

epvecd = logspace(-1,1,20);
epvecr = logspace(-1,1,20);
AN = [4,8,16,32,64,128];
BN = [4,8,16,32,64,128];
rbf = @(ep,r) exp(-(ep*r).^2);
rbfdx = @(ep,r,dx) -2*ep^2*dx.*exp(-(ep*r).^2);

alpha = 1;
FM = 21;
epvec = logspace(-2,1,FM);
%epvec = [];
epv = .8;%epv = 1e-5;
ind = 3;

warning off

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
couplewidth = 1; % How many points included in coupling
Fcstyle = 2; % 1 = Finite difference or 2 = meshfree

yf = @(x,t) exp(-t)*(1-0.5*(x(:,1).^2+x(:,2).^2));
ff = @(x,t) exp(-t)*(1+0.5*(x(:,1).^2+x(:,2).^2));

% yf = @(x,t) exp(-t)*sinh(x(:,1)+x(:,2));
% ff = @(x,t) -3*yf(x,t);

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

% The coupling style just matches values and derivatives
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
        end
        Frhs(FAifa) = zeros(size(FAifa));
        
        Fmat(FBifa,:) = zeros(size(Fmat(FBifa,:)));
        for k=1:length(FBifa)
            if couplewidth==1
                Fmat(FBifa(k),[FAcou(k),FAifa(k),FBcou(k),FBifa(k)]) = [[-1 1]/Adeltax [-1 1]/Bdeltax];
            elseif couplewidth==2
                Fmat(FBifa(k),[FAcou(2*k-[1,0]),FAifa(k),FBcou(2*k-[0,1]),FBifa(k)]) = [[-1 4 -3]/(2*Adeltax) [-1 4 -3]/(2*Bdeltax)];
            else
                error('couplewidth should be 1 or 2')
            end
        end
        Frhs(FBifa) = zeros(size(FBifa));
    elseif Fcstyle==2 % Meshfree coupling
        Fmat(FAifa,:) = zeros(size(Fmat(FAifa,:)));
        for k=1:length(FAifa)
            Fmat(FAifa(k),[FAifa(k),FBifa(k)]) = [1 -1];
        end
        Frhs(FAifa) = zeros(size(FAifa));
        Frhs(FBifa) = zeros(size(FBifa));
        
        AMFpts = Fx([FAcou,FAifa],:);
        ADM_int = DistanceMatrix(AMFpts,AMFpts);
        ADMdx_int = DistanceMatrix(Fx(FAifa,:),AMFpts);
        ADMdx_diff = DifferenceMatrix(Fx(FAifa,1),AMFpts(:,1));

        BMFpts = Fx([FBcou,FBifa],:);
        BDM_int = DistanceMatrix(BMFpts,BMFpts);
        BDMdx_int = DistanceMatrix(Fx(FBifa,:),BMFpts);
        BDMdx_diff = DifferenceMatrix(Fx(FBifa,1),BMFpts(:,1));
        
        Fval = zeros(FM,1);
        Fdir = zeros(FM,1);
        Fapx = zeros(FM,1);
        k = 1;
        for ep=epvec
            A_int = rbf(ep,ADM_int);
            A_x = rbfdx(ep,ADMdx_int,ADMdx_diff);
            B_int = rbf(ep,BDM_int);
            B_x = rbfdx(ep,BDMdx_int,BDMdx_diff);

            Fmat(FBifa,:) = zeros(size(Fmat(FBifa,:)));
            Fmat(FBifa,[FAcou,FAifa]) = A_x/A_int;
            Fmat(FBifa,[FBcou,FBifa]) = -B_x/B_int;
            Funew = Fmat\Frhs;
            Fdir(k) = errcompute(Funew,Fsol);
            
            ap = length(AMFpts);
%             Marr = rbfformMarr([0;0],[],2*ap);
%             Marr = [1 1 2 1 2 1 2 2  3 2 3 2 3 3 4 3;1 2 2 3 3 4 4 5  5 6 6 7 7 8 8 9];
%             Marr = [1 1 2 1 2 1 2 2  1 2 3 1 2 3 1 2;1 2 2 3 3 4 4 5  7 6 5 8 7 6 9 8];
%            Marr = rbfformMarr([couplewidth+1;ap]);
%            beta = (1+(2*ep/alpha)^2)^.25;
%            delta2 = .5*alpha^2*(beta^2-1);
%            Lvec = (alpha^2/(alpha^2+delta2+ep^2))* ...
%                      (ep^2/(alpha^2+delta2+ep^2)).^(sum(Marr)-2);
%            L1 = diag(Lvec(1:ap));
%            L2 = diag(Lvec(ap+1:end));
%            phiFull = rbfphi(Marr,AMFpts,ep,alpha);
%            phi1 = phiFull(:,1:ap);
%            phi2 = phiFull(:,ap+1:end);
%            phix = rbfphi(Marr,Fx(FAifa,:),ep,alpha,[1 0]);
%            phix1 = phix(:,1:ap);
%            phix2 = phix(:,ap+1:end);
%            A_psi = phiFull*[eye(ap);L2*(phi2'*pinv(phi1'))/L1];
%            A_psi_x = phix*[eye(ap);L2*(phi2'*pinv(phi1'))/L1];
%            
%            Fmat(FBifa,[FAcou,FAifa]) = A_psi_x*pinv(A_psi);
%            Fmat(FBifa,[FBcou,FBifa]) = A_psi_x*pinv(A_psi);
%            
%            Funew = Fmat\Frhs;
%            Fval(k) = errcompute(Funew,Fsol);
            
%             Marr = rbfformMarr([couplewidth+1;ap],[],floor(.7*ap));
            Marr = rbfformMarr([couplewidth+1;ap],[],ap);
            phiA = rbfphi(Marr,AMFpts,ep,alpha);
            phiAx = rbfphi(Marr,Fx(FBifa,:),ep,alpha,[1 0]);
            phiB = rbfphi(Marr,BMFpts,ep,alpha);
            phiBx = rbfphi(Marr,Fx(FBifa,:),ep,alpha,[1 0]);

            Fmat(FBifa,[FAcou,FAifa]) = phiAx/phiA;
            Fmat(FBifa,[FBcou,FBifa]) = -phiBx/phiB;
            
            Funew = Fmat\Frhs;
%            Funew = pinv(full(Fmat))*Frhs;
            Fapx(k) = errcompute(Funew,Fsol);

            fprintf('%d\t%g\t%g\t%g\n',k,ep,Fdir(k),Fapx(k))
            k = k+1;
        end
        if not(isempty(epvec)) % If we did a search for the best alpha
            [err,epi] = min(Fapx);
            epv = epvec(epi);
        end

        Marr = rbfformMarr([couplewidth+1;ap],[],ap);
        phiA = rbfphi(Marr,AMFpts,epv,alpha);
        phiAx = rbfphi(Marr,Fx(FBifa,:),epv,alpha,[1 0]);
        phiB = rbfphi(Marr,BMFpts,epv,alpha);
        phiBx = rbfphi(Marr,Fx(FBifa,:),epv,alpha,[1 0]);

        Fmat(FBifa,[FAcou,FAifa]) = phiAx/phiA;
        Fmat(FBifa,[FBcou,FBifa]) = -phiBx/phiB;
    end
    
    % Solve the system and compute the error on the grid
    Funew = Fmat\Frhs;
    Ferrt(tn) = errcompute(Funew,Fsol);
    
    % Unscatter the components back to A and B
    temp = Funew(Funshuffle);
    Aunew = temp(1:AMM);
    Bunew = temp(AMM+1:end);
end
warning on

% ex11
% This example is meshfree coupling between domains
rbfsetup
GAUSSQR_PARAMETERS.ERROR_STYLE = 2; % Absolute error

epvecd = logspace(-1,1,20);
epvecr = logspace(-1,1,20);
AN = [4,8,12,16,20];
BN = [4,8,12,16,20];
NN = 50;

rbf = @(ep,x) exp(-(ep*x).^2);
ep = 1;

dt = .1;
nsteps = 10;
stencil9 = [1/6 2/3 1/6  2/3 -10/3 2/3  1/6 2/3 1/6];
errt = zeros(1,nsteps);
couplebuffer = .1;

yf = @(x,t) exp(-t)*(1-.5*(x(:,1).^2+x(:,2).^2));
ff = @(x,t) exp(-t)*(1+.5*(x(:,1).^2+x(:,2).^2));

if not(exist('ind')) % Consider the first index (looping only)
  ind = 1;
end

Ald = [-1 0];Aud = [0 1];
Bld = [0 0];Bud = [1 1];
AM = AN(ind);
BM = BN(ind);

Axx = pick2Dpoints(Ald,Aud,NN*ones(2,1));
Bxx = pick2Dpoints(Bld,Bud,NN*ones(2,1));

Ax = pick2Dpoints(Ald,Aud,AN(ind)*ones(2,1));
Bx = pick2Dpoints(Bld,Bud,BN(ind)*ones(2,1));
AMM = size(Ax,1);
BMM = size(Bx,1);

%ishuffle = reshape(reshape(1:AMM,AM,AM)',AMM,1);
%Ax(:,:) = Ax(ishuffle,:)
%ishuffle = reshape(reshape(1:BMM,BM,BM)',BMM,1);
%Bx(:,:) = Bx(ishuffle,:)

Aifapts = find(abs(Ax(:,1)-0)<couplebuffer)';
Acoupts = setdiff(find(abs(Ax(:,1)-0)<=Adeltax*(1+couplebuffer)),Aifapts);
Aonlpts = setdiff(1:AMM,[Aifapts,Acoupts]);
Aintpts = find(1-any((Ax==1)+(Ax==0)+(Ax==-1),2))';
Aboupts = setdiff(1:AMM,Aintpts);

Ashuffle = [Aonlpts,Acoupts,Aifapts]
[temp,Aunshuffle] = sort(Ashuffle);
Ax = Ax(Ashuffle,:)
Astenoff9 = [-AM-1 -AM -AM+1  -1 0 1  AM-1 AM AM+1];
Asten = zeros(length(Aintpts),length(Astenoff9));
for k=1:length(Aintpts)
  Asten(k,:) = Aunshuffle(Astenoff9+Aintpts(k));
end
Adeltax = 1/(AM-1);

Aifapts = find(abs(Ax(:,1)-0)<couplebuffer)';
Acoupts = setdiff(find(abs(Ax(:,1)-0)<=Adeltax*(1+couplebuffer)),Aifapts);
Aonlpts = setdiff(1:AMM,[Aifapts,Acoupts]);
Aintpts = find(1-any((Ax==1)+(Ax==0)+(Ax==-1),2))';
Aboupts = setdiff(1:AMM,Aintpts);

[temp,order] = sort(Asten(:,ceil(length(Asten)/2)));
Asten = Asten(order,:);

pause

Bintpts = find(1-any((Bx==1)+(Bx==0)+(Bx==-1),2))';
Bboupts = setdiff(1:BMM,intpts);
Bifapts = find(abs(Bx(:,1)-0)<couplebuffer);
Bcoupts = setdiff(find(abs(Bx(:,1)-0)<=Bdeltax*(1+couplebuffer)),Bifapts);
Bstenoff9 = [-BM-1 -BM -BM+1  -1 0 1  BM-1 BM BM+1];
Bdeltax = 1/(BM-1);

Aunew = yf(Ax,0);
Bunew = yf(Bx,0);

for tn=1:nsteps
  t = tn*dt;

  Amat = 1/dt*eye(AMM);
  Arhs = zeros(AMM,1);
  Arf = ff(Ax,t);
  Aubc = yf(Ax,t);
  Auold = Aunew;

  m = 1;
  for k=Aintpts
      Amat(k,Asten(m,:)) = Amat(k,Asten(m,:)) - stencil9/(Adeltax^2);
%      Amat(k,Astenoff9+k) = Amat(k,Astenoff9+k) - stencil9/(Adeltax^2);
      Arhs(k) = Arf(k) + 1/dt*Auold(k);
      m = m + 1;
  end
  for k=Aboupts
      Amat(k,k) = 1.0;
      Arhs(k) = Aubc(k);
  end

  Aunew = Amat\Arhs;
  Aerrt(tn) = errcompute(Aunew,Aubc);


  Bmat = 1/dt*eye(BMM);
  Brhs = zeros(BMM,1);
  Brf = ff(Bx,t);
  Bubc = yf(Bx,t);
  Buold = Bunew;

  for k=Bintpts
      Bmat(k,Bstenoff9+k) = Bmat(k,Bstenoff9+k) - stencil9/(Bdeltax^2);
      Brhs(k) = Brf(k) + 1/dt*Buold(k);
  end
  for k=Bboupts
      Bmat(k,k) = 1.0;
      Brhs(k) = Bubc(k);
  end

  Bunew = Bmat\Brhs;
  Berrt(tn) = errcompute(Bunew,Bubc);

  Fmat = [Amat,zeros(BMM,AMM);zeros(AMM,BMM),Bmat];
end

% ex11b
% Testing the validity of derivatives for psi and phi functions at the
% interface manifold for a coupling problem
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2; % Absolute error

AN = 13;
BN = AN;
rbf = @(ep,r) exp(-(ep*r).^2);
rbfdx = @(ep,r,dx) -2*ep^2*dx.*exp(-(ep*r).^2);

alpha = 3;
epN = 25;
epvec = logspace(-5,.5,epN);
Dir_x = zeros(1,epN);
QR_x = zeros(1,epN);
QRr_x = zeros(1,epN);

yf = @(x) 1-0.5*(x(:,1).^2+x(:,2).^2);
yfdx = @(x) -2*x(:,1);
yf = @(x) sin(x(:,1).*x(:,2));
yfdx = @(x) x(:,2).*cos(x(:,1).*x(:,2));
% yf = @(x) exp(x(:,1)-x(:,2));
% yfdx = @(x) exp(x(:,1)-x(:,2));

Ald = [-1 0];Aud = [0 1];
Bld = [0 0];Bud = [1 1];
couplebuffer = 0.1/AN; % accounts for fuzzy math
couplewidth = 3; % How many points included in coupling

% All the points in the problem, even the non coupling ones
Ax = pick2Dpoints(Ald,Aud,AN*ones(2,1));
Adeltax = 1/(AN-1);
Bx = pick2Dpoints(Bld,Bud,BN*ones(2,1));
Bdeltax = 1/(BN-1);

Aifapts = find(abs(Ax(:,1)-0)<couplebuffer)';
Acoupts = setdiff(find(abs(Ax(:,1)-0)<=Adeltax*(couplewidth+couplebuffer)),Aifapts);

Bifapts = find(abs(Bx(:,1)-0)<couplebuffer)';
Bcoupts = setdiff(find(abs(Bx(:,1)-0)<=Bdeltax*(couplewidth+couplebuffer)),Bifapts);

Fx = [Ax([Acoupts,Aifapts],:);Bx([Bcoupts,Bifapts],:)];

FAind = 1:length([Acoupts,Aifapts]);
FBind = FAind(end)+1:FAind(end)+length([Bcoupts,Bifapts]);
FAcou = 1:length(Acoupts);
FAifa = FAcou(end)+1:FAcou(end)+length(Aifapts);
FBcou = FAifa(end)+1:FAifa(end)+length(Bcoupts);
FBifa = FBcou(end)+1:FBcou(end)+length(Bifapts);

yFx = yf(Fx(FAind,:));
true_x = yfdx(Fx(FBifa,:));

% NOTE: Need to fix this so that model A gets evaluated at B coupling
% locations
AMFpts = Fx(FAind,:);
ADM_int = DistanceMatrix(AMFpts,AMFpts);
ADMdx_int = DistanceMatrix(Fx(FBifa,:),AMFpts);
ADMdx_diff = DifferenceMatrix(Fx(FBifa,1),AMFpts(:,1));

BMFpts = Fx(FBind,:);
BDM_int = DistanceMatrix(BMFpts,BMFpts);
BDMdx_int = DistanceMatrix(Fx(FBifa,:),BMFpts);
BDMdx_diff = DifferenceMatrix(Fx(FBifa,1),BMFpts(:,1));

% Dmat is the part of the matrix which dealing with the exchanging of
% derivatives between models.  The assumption here is that model A is
% calculating its derivatives on model B's points.  That's why the size of
% the matrix conforms to model B
Dmat = zeros(length(FBifa),length(FAind));

warning off
k = 1;
for ep=epvec
    A_int = rbf(ep,ADM_int);
    A_x = rbfdx(ep,ADMdx_int,ADMdx_diff);
%     B_int = rbf(ep,BDM_int);
%     B_x = rbfdx(ep,BDMdx_int,BDMdx_diff);

%     Dmat(:,[FAcou,FAifa]) = A_x/A_int;
%     Dmat(:,[FBcou,FBifa]) = -B_x/B_int;
    Dmat = A_x/A_int;
    PS_dir = Dmat*yFx;
    Dir_x(k) = errcompute(PS_dir,true_x);

    ap = length(FAind);
    bp = length(FBind);
%     Marr = gqr_formMarr([0;0],[],ap+bp);
%     Marr = [1 1 2 1 2 1 2 2  1 2 3 1 2 3 1 2;1 2 2 3 3 4 4 5  7 6 5 8 7 6 9 8];
    Marr = gqr_formMarr([couplewidth+1;ap]);
    beta = (1+(2*ep/alpha)^2)^.25;
    delta2 = .5*alpha^2*(beta^2-1);
    Lvec = (alpha^2/(alpha^2+delta2+ep^2))* ...
              (ep^2/(alpha^2+delta2+ep^2)).^(sum(Marr)-2);
    L1 = diag(Lvec(1:ap));
    L2 = diag(Lvec(ap+1:end));
    
    phiA = gqr_phi(Marr,AMFpts,ep,alpha);
    phiA1 = phiA(:,1:ap);
    phiA2 = phiA(:,ap+1:end);
    phiAx = gqr_phi(Marr,Fx(FBifa,:),ep,alpha,[1 0]);
    phiAx1 = phiAx(:,1:ap);
    phiAx2 = phiAx(:,ap+1:end);
    Abar = [eye(ap);L2*(phiA2'*pinv(phiA1'))/L1];
    A_psi = phiA*Abar;
    A_psi_x = phiAx*Abar;
    
%     phiB = gqr_phi(Marr,BMFpts,ep,alpha);
%     phiB1 = phiB(:,1:bp);
%     phiB2 = phiB(:,bp+1:end);
%     phiBx = gqr_phi(Marr,Fx(FBifa,:),ep,alpha,[1 0]);
%     phiBx1 = phiBx(:,1:bp);
%     phiBx2 = phiBx(:,bp+1:end);
%     Bbar = [eye(ap);L2*(phiB2'*pinv(phiB1'))/L1];
%     B_psi = phiB*Bbar;
%     B_psi_x = phiBx*Bbar;

%     Dmat(:,[FAcou,FAifa]) = A_psi_x/A_psi;
%     Dmat(:,[FBcou,FBifa]) = -B_psi_x/B_psi;
    Dmat = A_psi_x/A_psi;
    PS_QR = Dmat*yFx;
    QR_x(k) = errcompute(PS_QR,true_x);
    
    Marr = gqr_formMarr([couplewidth+1;ap],[],floor(.7*ap));
%     Marr = gqr_formMarr([0;0],[],floor(.7*ap));
    phiA = gqr_phi(Marr,AMFpts,ep,alpha);
    phiAx = gqr_phi(Marr,Fx(FBifa,:),ep,alpha,[1 0]);
%     phiB = gqr_phi(Marr,BMFpts,ep,alpha);
%     phiBx = gqr_phi(Marr,Fx(FBifa,:),ep,alpha,[1 0]);

%     Dmat(:,[FAcou,FAifa]) = phiAx/phiA;
%     Dmat(:,[FBcou,FBifa]) = -phiBx/phiB;
    Dmat = phiAx/phiA;
    PS_QRr = Dmat*yFx;
    QRr_x(k) = errcompute(PS_QRr,true_x);
    
    k = k+1;
end
warning on

loglog(epvec,[Dir_x;QR_x;QRr_x])
% ex19c
% This example uses the multiquadrics on the sphere

gqr_downloaddata('sphereMDpts_data.mat')
load sphereMDpts_data

yf = @(x) (1 + 10*x(:,1).*x(:,2).*x(:,3) + x(:,1).^8 + exp(2*x(:,2).^3) + exp(2*x(:,3).^2))/14 - 3*cos(7*x(:,1)-2*x(:,3)) + 2*exp(-2*x(:,2).^2);

x = sphereMDpts{3};
N = size(x,1);
y = yf(x);
DP = x*x';

NN = 2000;
[t,testptstr] = pick2Dpoints([-pi -pi/2],[pi pi/2],sqrt(NN),'halt');
xeval = zeros(NN,3);
[xeval(:,1),xeval(:,2),xeval(:,3)] = sph2cart(t(:,1),t(:,2),1);
yy = yf(xeval);
DPeval = xeval*x';

zbf = @(g,dp) 1./sqrt(1+g^2-2*g*dp);

gammavec = logspace(-2,-.1,30);

Nlevels = ceil(sqrt(N)) + 3;
Neig = Nlevels^2;
lrange = 0:Nlevels-1;
larr = cell2mat(arrayfun(@(l)l*ones(1,2*l+1),lrange,'UniformOutput',0));
marr = cell2mat(arrayfun(@(l)1:2*l+1,lrange,'UniformOutput',0));
Phi = cell2mat(arrayfun(@(l,m)sphHarm(l,m,x),larr,marr,'UniformOutput',0));
Phieval = cell2mat(arrayfun(@(l,m)sphHarm(l,m,xeval),larr,marr,'UniformOutput',0));

% Break Phi and Phieval into the blocks for the HSSVD
Phi1 = Phi(:,1:N);
Phi2 = Phi(:,N+1:Neig);
Phieval1 = Phieval(:,1:N);
Phieval2 = Phieval(:,N+1:Neig);
Phi2TinvPhi1 = Phi2'/Phi1';

% Perform the interpolation with the range of gamma and record the values
k = 1;
errvec = zeros(size(gammavec));
errvechs = zeros(size(gammavec));
mlevec = zeros(size(gammavec));
mlevechs = zeros(size(gammavec));
for g=gammavec
    % Perform the standard interpolation
    K = zbf(g,DP);
    Keval = zbf(g,DPeval);
    warning('off','MATLAB:nearlySingularMatrix')
    c = K\y;
    warning('on','MATLAB:nearlySingularMatrix')
    yeval = Keval*c;
    
    % Eigs are the eigenvalues of the kernel, here stored as a vector since
    % storing them as a diagonal matrix seems unnecessary.
    % We create LamFull which is the matrix that has the effect of both
    % multiplying on the left by Lam_2 and on the right by inv(Lam_1)
    lamvec = 4*pi*g.^larr./(2*larr+1);
    lamvec1 = lamvec(1:N);
    lamvec2 = lamvec(N+1:Neig);
    LamFull = bsxfun(@rdivide,lamvec2',lamvec1);
    CbarT = LamFull.*Phi2TinvPhi1;
    Psi = Phi1 + Phi2*CbarT;
    b = Psi\y;
    Psieval = Phieval1 + Phieval2*CbarT;
    yevalhs = Psieval*b;
    
    logdetPsi = sum(log(svd(Psi)));
    logdetPhi = sum(log(svd(Phi1)));
    logdetLam = sum(log(lamvec1));
    logdetK = logdetPsi + logdetPhi + logdetLam;
    boundvec = b'*(b./lamvec1');
    L2P2P1L1invb = (lamvec2'.^-.5).*(CbarT*b);
    correctionvec = L2P2P1L1invb'*L2P2P1L1invb;
    mahaldist = boundvec + correctionvec;
    
    errvechs(k) = errcompute(yevalhs,yy);
    errvec(k) = errcompute(yeval,yy);
    mlevechs(k) = N*log(mahaldist) + logdetK;
    mlevec(k) = N*log(abs(c'*y)) + sum(log(svd(K)));
    k = k + 1;
end

% Plot the results
h_joint = figure;

yyaxis left
loglog(gammavec,errvechs,'linewidth',3)
hold on
loglog(gammavec,errvec,'--','linewidth',2)
hold off
xlabel('\gamma','fontsize',14)
ylabel('pointwise error','fontsize',14)
ax = gca;
ax.YTick = [1e-5, 1, 1e4];

yyaxis right
semilogx(gammavec,mlevechs,'linewidth',3)
hold on
semilogx(gammavec,mlevec,'--','linewidth',2)
hold off
ylabel('likelihood','fontsize',14)
ax = gca;
ax.YTick = [-2000, 0, 4000];
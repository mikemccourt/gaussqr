clear all

delta_ep = 0.1;
epsilon = 0.01:delta_ep:50;
% epsilon = logspace(-2,2,50);

N=size(epsilon,2);
rel_err_phi = zeros(N,1);
COND = zeros(N,1);

% Medium data
R = 0.1;                      % Sphere radius [m]
sig = 0.2;                    % Electric conductivity [S/m]
% Sources data
dipmom = 2.7e-12.*[1, 0, 0];  % Dipole moment [Am]
srcpnts = [0, 0, 0.6*R];      % Dipole position [m]
% Magnetic induction field observation points
% obspnts = ; ...
% Parameters for numerical computation
Npnts = 700;                % Number of desired interior points

rbfset = { 'Gaussian' 'IMQ' 'MQ' 'LinearMatern' ...
           'Wendland_C2' 'Wendland_C4' 'Wendland_C6' };
% rbfset = { 'Gaussian' };

% Collocation points
[POINTS, NORMALS] = BallGeometry(R, Npnts, 'kansa');
intdata = POINTS.int1;
bdydata = POINTS.bdy11;

Npnts_int = length(intdata);
Npnts_bdy = length(bdydata);

% Boundary centers (inside the domain)
bdyctrs = bdydata;

% Centers
ctrs = [intdata; bdydata];

% Evaluation points
evalpnts = SphereSurfGoldPoints(1000, R);
neval = size(evalpnts,1);

% Compute evaluation matrix EM
DM_eval = DistanceMatrix(evalpnts, ctrs);

% Compute blocks for collocation matrix
% Interior points
DM_intdata = DistanceMatrix(intdata,ctrs);

% Boundary points
DM_bdydata = DistanceMatrix(bdydata,ctrs);
dx_bdydata = DifferenceMatrix(bdydata(:,1),ctrs(:,1));
dy_bdydata = DifferenceMatrix(bdydata(:,2),ctrs(:,2));
dz_bdydata = DifferenceMatrix(bdydata(:,3),ctrs(:,3));

% Potential at evalpnts in the unbound domain case
phi_F = phiF_dip(evalpnts,srcpnts,dipmom,sig);

phi_an = HomSpherePotential(R, sig, srcpnts, dipmom, evalpnts);

% Compute known-terms vector (a.k.a. righthand side vector)
gradphi_F = gradphiF_dip(bdydata, srcpnts, dipmom, sig);% Gradient of the 
                                                    % potential at boundary
                                                    % in the unbound case
rhs = [ zeros(size(intdata,1),1); -sum(NORMALS.n11.*gradphi_F,2) ];


for j = 1:length(rbfset)
[rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF(rbfset{j});

for i=1:N
    [rel_err_phi(i,1), COND(i) ] = ...
        epskansa_dip(epsilon(i), NORMALS, rhs, DM_eval, DM_intdata, DM_bdydata, ...
                 dx_bdydata, dy_bdydata, dz_bdydata, phi_F, phi_an, ...
                 rbf, dxrbf, dyrbf, dzrbf, Lrbf);
    fprintf('%1s RBF, epsilon %4u out of %4u\n', rbfset{j}, i, N)
end

h = figure('Name', strcat(rbfset{j}, ' RBF, ', ...
                          num2str(Npnts_int), ' interior points, ', ...
                          num2str(Npnts_bdy), ' boundary points'));
subplot(2,1,1)
loglog(epsilon, rel_err_phi);
xlabel('\epsilon');
ylabel('Relative Error');
subplot(2,1,2)
loglog(epsilon, COND);
xlabel('\epsilon');
ylabel('cond(CM)');

filename = strcat(rbfset{j},'_',num2str(Npnts_int),'int_',num2str(Npnts_bdy),'bdy');
saveas(h, filename,'fig')
save(filename, 'rel_err_phi', 'COND', 'epsilon')

end
clear all

a=0;
delta_ep = 0.05;
epsilon = [0.05:delta_ep:20];

N=size(epsilon,2)
rms_err_phi = zeros(N,1);
rel_err_phi = zeros(N,1);
max_err_phi = zeros(N,1);
COND = zeros(N,1);


% Medium data
R = 0.1;                      % Sphere radius [m]
sig = 0.2;                    % Electric conductivity [S/m]
% Sources data
dipmom = 2.7e-12.*[1, 0, 0];         % Dipole moment [Am]
srcpnts = [0, 0, 0.6*R];  % Dipole position [m]
% Parameters for numerical computation
Npnts_surf = 300;      % Number of desired points on sphere's surface 300

% Collocation points
[bdydata, intdata, dist] = SphereRegGoldPoints(Npnts_surf, R);
% Boundary centers outside the domain
bdyctrs = SphereSurfGoldPoints(Npnts_surf, R+a*dist);
% Centers
ctrs = [intdata; bdyctrs];
% Evaluation points
evalpnts = SphereRegUnifPoints(dist/3, R); % Dist/3 controls the distance 
                                           % between evaluation points
neval = size(evalpnts,1);

% Compute known-terms vector (a.k.a. righthand side vector)
gradphi_F = gradphiF(bdydata, srcpnts, dipmom, sig);% Gradient of the 
                                                    % potential at boundary
                                                    % in the unbound case
NV = bdydata/R; % Unit vectors normal to sphere surface
rhs = [ zeros(size(intdata,1),1); -sum(NV.*gradphi_F,2) ];

% Potential at evalpnts in the unbound domain case
phi_F = phiF(evalpnts,srcpnts,dipmom,sig);
%  Analytic solution for the potential
phi_an = HomSpherePotential(R, sig, srcpnts, dipmom, evalpnts);

DM_eval = DistanceMatrix(evalpnts, ctrs);
DM_intdata = DistanceMatrix(intdata,ctrs);
DM_bdydata = DistanceMatrix(bdydata,ctrs);
dx_bdydata = DifferenceMatrix(bdydata(:,1),ctrs(:,1));
dy_bdydata = DifferenceMatrix(bdydata(:,2),ctrs(:,2));
dz_bdydata = DifferenceMatrix(bdydata(:,3),ctrs(:,3));

% %% IMQ, bdycenters on the boundary, 300 pnts on surface
% [rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF('imq');
% for i=1:N
%     [rms_err_phi(i,1), rel_err_phi(i,1), max_err_phi(i,1), COND(i) ] = epskansa(epsilon(i), NV, rhs, DM_eval, DM_intdata, DM_bdydata, dx_bdydata, dy_bdydata, dz_bdydata, phi_F, phi_an, rbf, dxrbf, dyrbf, dzrbf, Lrbf, neval);
%     i
% end
% 
% figure('Name','IMQ, 300 pnts on surface');
% 
% subplot(2,2,1)
% plot(epsilon, rel_err_phi);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('Relative Error');
% 
% subplot(2,2,2)
% plot(epsilon, rms_err_phi);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('RMS Error');
% 
% subplot(2,2,3)
% plot(epsilon, max_err_phi);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('Max Error');
% 
% subplot(2,2,4)
% plot(epsilon, COND);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('cond(CM)');
% 
% save('IMQ_300pnts')
% 
% 
% %% Gaussian, bdycenters on the boundary, 300 pnts on surface
% 
% [rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF('gaussian');
% 
% for i=1:N
%     [rms_err_phi(i,1), rel_err_phi(i,1), max_err_phi(i,1), COND(i) ] = epskansa(epsilon(i), NV, rhs, DM_eval, DM_intdata, DM_bdydata, dx_bdydata, dy_bdydata, dz_bdydata, phi_F, phi_an, rbf, dxrbf, dyrbf, dzrbf, Lrbf, neval);
%     i
% end
% 
% figure('Name','Gaussian, 300 pnts on surface');
% 
% subplot(2,2,1)
% plot(epsilon, rel_err_phi);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('Relative Error');
% 
% subplot(2,2,2)
% plot(epsilon, rms_err_phi);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('RMS Error');
% 
% subplot(2,2,3)
% plot(epsilon, max_err_phi);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('Max Error');
% 
% subplot(2,2,4)
% plot(epsilon, COND);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('cond(CM)');
% 
% save('Gaussian_300pnts')
% 
% %% Matern linear, bdycenters on the boundary, 300 pnts on surface
% 
% [rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF('linearmatern');
% 
% for i=1:N
%     [rms_err_phi(i,1), rel_err_phi(i,1), max_err_phi(i,1), COND(i) ] = epskansa(epsilon(i), NV, rhs, DM_eval, DM_intdata, DM_bdydata, dx_bdydata, dy_bdydata, dz_bdydata, phi_F, phi_an, rbf, dxrbf, dyrbf, dzrbf, Lrbf, neval);
%     i
% end
% 
% figure('Name','Matern (linear), 300 pnts on surface');
% 
% subplot(2,2,1)
% plot(epsilon, rel_err_phi);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('Relative Error');
% 
% subplot(2,2,2)
% plot(epsilon, rms_err_phi);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('RMS Error');
% 
% subplot(2,2,3)
% plot(epsilon, max_err_phi);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('Max Error');
% 
% subplot(2,2,4)
% plot(epsilon, COND);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('cond(CM)');
% 
% save('maternlin_300pnts')
% 
% 
% %% MQ, bdycenters on the boundary, 300 pnts on surface
% 
% [rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF('mq');
% 
% for i=1:N
%     [rms_err_phi(i,1), rel_err_phi(i,1), max_err_phi(i,1), COND(i) ] = epskansa(epsilon(i), NV, rhs, DM_eval, DM_intdata, DM_bdydata, dx_bdydata, dy_bdydata, dz_bdydata, phi_F, phi_an, rbf, dxrbf, dyrbf, dzrbf, Lrbf, neval);
%     i
% end
% 
% figure('Name','MQ, 300 pnts on surface');
% 
% subplot(2,2,1)
% plot(epsilon, rel_err_phi);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('Relative Error');
% 
% subplot(2,2,2)
% plot(epsilon, rms_err_phi);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('RMS Error');
% 
% subplot(2,2,3)
% plot(epsilon, max_err_phi);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('Max Error');
% 
% subplot(2,2,4)
% plot(epsilon, COND);
% set(gca,'Yscale','log');
% xlabel('\epsilon');
% ylabel('cond(CM)');
% 
% save('MQ_300pnts')

%% Wendland C2, bdycenters on the boundary, 300 pnts on surface

[rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF('wendland_c2');

for i=1:N
    [rms_err_phi(i,1), rel_err_phi(i,1), max_err_phi(i,1), COND(i) ] = epskansa(epsilon(i), NV, rhs, DM_eval, DM_intdata, DM_bdydata, dx_bdydata, dy_bdydata, dz_bdydata, phi_F, phi_an, rbf, dxrbf, dyrbf, dzrbf, Lrbf, neval);
    i
end

figure('Name','Wendland C2, 300 pnts on surface');

subplot(2,2,1)
plot(epsilon, rel_err_phi);
set(gca,'Yscale','log');
xlabel('\epsilon');
ylabel('Relative Error');

subplot(2,2,2)
plot(epsilon, rms_err_phi);
set(gca,'Yscale','log');
xlabel('\epsilon');
ylabel('RMS Error');

subplot(2,2,3)
plot(epsilon, max_err_phi);
set(gca,'Yscale','log');
xlabel('\epsilon');
ylabel('Max Error');

subplot(2,2,4)
plot(epsilon, COND);
set(gca,'Yscale','log');
xlabel('\epsilon');
ylabel('cond(CM)');

save('WC2_300pnts')

%% Wendland C4, bdycenters on the boundary, 300 pnts on surface

[rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF('wendland_c4');

for i=1:N
    [rms_err_phi(i,1), rel_err_phi(i,1), max_err_phi(i,1), COND(i) ] = epskansa(epsilon(i), NV, rhs, DM_eval, DM_intdata, DM_bdydata, dx_bdydata, dy_bdydata, dz_bdydata, phi_F, phi_an, rbf, dxrbf, dyrbf, dzrbf, Lrbf, neval);
    i
end

figure('Name','Wendland C4, 300 pnts on surface');

subplot(2,2,1)
plot(epsilon, rel_err_phi);
set(gca,'Yscale','log');
xlabel('\epsilon');
ylabel('Relative Error');

subplot(2,2,2)
plot(epsilon, rms_err_phi);
set(gca,'Yscale','log');
xlabel('\epsilon');
ylabel('RMS Error');

subplot(2,2,3)
plot(epsilon, max_err_phi);
set(gca,'Yscale','log');
xlabel('\epsilon');
ylabel('Max Error');

subplot(2,2,4)
plot(epsilon, COND);
set(gca,'Yscale','log');
xlabel('\epsilon');
ylabel('cond(CM)');

save('WC4_300pnts')


%% Wendland C6, bdycenters on the boundary, 300 pnts on surface

[rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF('wendland_c6');

for i=1:N
    [rms_err_phi(i,1), rel_err_phi(i,1), max_err_phi(i,1), COND(i) ] = epskansa(epsilon(i), NV, rhs, DM_eval, DM_intdata, DM_bdydata, dx_bdydata, dy_bdydata, dz_bdydata, phi_F, phi_an, rbf, dxrbf, dyrbf, dzrbf, Lrbf, neval);
    i
end

figure('Name','Wendland C6, 300 pnts on surface');

subplot(2,2,1)
plot(epsilon, rel_err_phi);
set(gca,'Yscale','log');
xlabel('\epsilon');
ylabel('Relative Error');

subplot(2,2,2)
plot(epsilon, rms_err_phi);
set(gca,'Yscale','log');
xlabel('\epsilon');
ylabel('RMS Error');

subplot(2,2,3)
plot(epsilon, max_err_phi);
set(gca,'Yscale','log');
xlabel('\epsilon');
ylabel('Max Error');

subplot(2,2,4)
plot(epsilon, COND);
set(gca,'Yscale','log');
xlabel('\epsilon');
ylabel('cond(CM)');

save('WC6_300pnts')
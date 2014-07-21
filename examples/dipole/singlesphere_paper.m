%% singlesphere_paper.m
% This script runs the simulation we need for the paper (homogeneous
% sphere case): it is based on singlesphere_N.m, so look at the help of
% that file for further information.
% Here we investigate convergence, conditioning and cost of our MFS solver
% (with different mfs_frac) and make comparisons with Kansa and BEM.
%

% Choose the solvers you want to use
BEM = 0;
MFS = 1;
kansa = 1;

R = 0.1;
sig = 0.02;
dipmom = [1, 0, 0];
srcpnts = [0, 0, 0.6*R];

radbasfun = 'imq';
ep = 10;

int_point_dist = 'halton';
bdy_point_dist = 'spiral';

mfs_frac = [1.0 0.8 0.4 0.2];
mfs_sphere = 1.3;

reference = [0, 0, -R];

Nvec = 350:500:5350;
N_eval = 1000;

iter_out = 1;
plot_sol = 1;
sol_err_style = 1;

% N_plotsol = length(Nvec);
N_plotsol = 3;
mfs_frac_plotsol = 1;

resultsfilename = 'homogeneous_sphere_tests';

%%%%%%%%%%%%%%%%%%%%%
% Basic setup stuff independent of this problem

% Consider the standard GQR parameters for the errcompute function
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
GAUSSQR_PARAMETERS.ERROR_STYLE = 3;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

% Set random number generator to constant
% This is used in choosing which BC points are Dirichlet in the mixed case
rng(0);


%%%%%%%%%%%%%%%%%%%%%
% This is the start of the solver

% Determine the evaluation points (all on the boundary)
evalpnts = SphereSurfGoldPoints(N_eval-1, R);
evalpnts = [reference;evalpnts];

% Potential at evalpnts in the unbound domain case
% This is the analytic component of the computed solution
tic;
phi_F = phiF_dip(evalpnts,srcpnts,dipmom,sig);
tt = toc;

% Analytic solution for the potential
phi_an = HomSpherePotential(R, sig, srcpnts, dipmom, evalpnts);

% Compute the difference of the solution with a reference
% point, arbitrarily chosen as evalpnts(1)
phi_true = phi_an - phi_an(1);

%% BEM solutions
if BEM
    N_elements = zeros(length(Nvec),1); 
    phi_comp_BEM = zeros(N_eval,length(Nvec));
    for k=1:length(Nvec)
        if iter_out
            fprintf('k=%d\n',k)
        end
        [vol, N_elements(k,:)] = prepare_multisphere_mesh(Nvec(k),R);
        [phi_comp_BEM(:,k),~,~,elapsed_t_BEM(k)] = BEM_potential(vol,sig,evalpnts,srcpnts,dipmom);
        errvec_BEM(k) = errcompute(phi_comp_BEM(:,k),phi_true);
        Nvec_true_BEM(k) = sum(N_elements(k,:));
        if iter_out
            fprintf('\terr = %g\n\tN = %d\n',errvec_BEM(k),Nvec_true_BEM(k));
        end
    end
end

%% MFS solutions
if MFS
    % RBF definition and derivatives
    [rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF('fundamental_3d');
    
    phi_comp_MFS = zeros(N_eval,length(Nvec),length(mfs_frac));
    for i=1:length(mfs_frac)
        
        for k=1:size(N_elements,1)
            
            Npnts = N_elements(k,:);
            
            if iter_out
                fprintf('k=%d\n',k)
            end
            
            % Determine collocation points
            [POINTS, NORMALS] = BallGeometry(R,Npnts,'mfs',int_point_dist,bdy_point_dist);
            bdydata = POINTS.bdy11;
            normvecs = NORMALS.n11;
            N_bdy = size(bdydata,1);
            N_tot = N_bdy;  % Total number of collocation points
            
            % Compose a vector of all the RBF centers
            % In the MFS setting, these are chosen in a sphere around the ball
            N_ctrs = floor(mfs_frac(i)*sum(Npnts));
            ctrs = SphereSurfGoldPoints(N_ctrs, mfs_sphere*R);
            
            tic;
            
            % Compute the evaluation matrix
            DM_eval = DistanceMatrix(evalpnts, ctrs);
            EM = rbf(ep, DM_eval);
            
            % Compute the collocation block for the boundary conditions
            % This also computes the RHS for the problem
            % First we consider the Neumann BC component
            DM_bdydata = DistanceMatrix(bdydata,ctrs);
            
            % Find all the necessary difference matrices
            dx_bdydata = DifferenceMatrix(bdydata(:,1),ctrs(:,1));
            dy_bdydata = DifferenceMatrix(bdydata(:,2),ctrs(:,2));
            dz_bdydata = DifferenceMatrix(bdydata(:,3),ctrs(:,3));
            
            % Compute normal derivative collocation matrix for boundary
            A = repmat(normvecs(:,1),1,N_ctrs).*dxrbf(ep,DM_bdydata,dx_bdydata);
            B = repmat(normvecs(:,2),1,N_ctrs).*dyrbf(ep,DM_bdydata,dy_bdydata);
            C = repmat(normvecs(:,3),1,N_ctrs).*dzrbf(ep,DM_bdydata,dz_bdydata);
            CM = A + B + C;
            
            % Compute known-terms vector (a.k.a. righthand side vector)
            % This requires the gradient of the unbounded potential at boundary
            gradphi_F_neu = gradphiF_dip(bdydata, srcpnts, dipmom, sig);
            rhs = -sum(normvecs.*gradphi_F_neu,2);
            
            %         [coefs,recip_cond] = linsolve(CM,rhs);
            [U,S,V] = svd(CM,0);
            sing_val = diag(S);
            inv_sing_val_TSVD = 1./sing_val;
            indices = find(sing_val < 1e-6*max(sing_val));
            inv_sing_val_TSVD(indices) = 0;
            % inv_sing_val_TSVD(floor(0.2*length(inv_sing_val)):end) = 0;
            cbeta = U'*rhs;
            coefs = V*(inv_sing_val_TSVD.*cbeta(1:length(sing_val)));
            
            % Potential at evalpnts in the source free case
            phi0 = EM * coefs;
            % Potential at evalpnts (superposition of effects)
            phi = phi0 + phi_F;
            
            % Compute the difference of the solution with a
            % reference point, arbitrarily chosen as evalpnts(1)
            phi_comp_MFS(:,k,i) = phi - phi(1);
            
            elapsed_t_MFS(k,i) = tt + toc;
            
            % Compute the total errors
            errvec_MFS(k,i) = errcompute(phi_comp_MFS(:,k,i),phi_true);
            
            % Store the condition of the system
            %         [U,S,~] = svd(CM);
            %         sing_val = diag(S);
            %         cbeta = U'*rhs;
            condvec_MFS(k,i) = norm(rhs)/min(sing_val)/...
                norm(cbeta(1:length(sing_val))./sing_val);
            condstr = 'Effective condition number';
            
            Nvec_true_MFS(k) = N_tot;
            
            if iter_out
                fprintf('\terr = %g\n\tcond = %g\n\tN = %d\n',errvec_MFS(k,i),condvec_MFS(k,i),Nvec_true_MFS(k));
            end
        end
        
    end
end
%% Kansa solutions
if kansa
    % RBF definition and derivatives
    [rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF(radbasfun);
    
    phi_comp_kansa = zeros(N_eval,length(Nvec));
    for k=1:length(Nvec)
        
        if iter_out
            fprintf('k=%d\n',k)
        end
        
        % Determine collocation points
        [POINTS, NORMALS] = BallGeometry(R,N_elements(k,:),'kansa',int_point_dist,bdy_point_dist);
        intdata = POINTS.int1;
        bdydata = POINTS.bdy11;
        normvecs = NORMALS.n11;
        N_int = size(intdata,1);
        N_bdy = size(bdydata,1);
        N_tot = N_int + N_bdy;  % Total number of collocation points
        
        if N_tot < max(Nvec)*1.15
            % Compose a vector of all the RBF centers
            ctrs = [intdata; bdydata];
            N_ctrs = N_tot;
            
            tic;
            % Compute the collocation block for the interior
            DM_intdata = DistanceMatrix(intdata,ctrs);
            LCM = Lrbf(ep,DM_intdata);
            rhs_int = zeros(N_int,1);
            
            % Compute the evaluation matrix
            DM_eval = DistanceMatrix(evalpnts, ctrs);
            EM = rbf(ep, DM_eval);
            
            % Compute the collocation block for the boundary conditions
            % This also computes the RHS for the problem
            % First we consider the Neumann BC component
            DM_bdydata = DistanceMatrix(bdydata,ctrs);
            
            % Find all the necessary difference matrices
            dx_bdydata = DifferenceMatrix(bdydata(:,1),ctrs(:,1));
            dy_bdydata = DifferenceMatrix(bdydata(:,2),ctrs(:,2));
            dz_bdydata = DifferenceMatrix(bdydata(:,3),ctrs(:,3));
            
            % Compute normal derivative collocation matrix for boundary
            A = repmat(normvecs(:,1),1,N_ctrs).*dxrbf(ep,DM_bdydata,dx_bdydata);
            B = repmat(normvecs(:,2),1,N_ctrs).*dyrbf(ep,DM_bdydata,dy_bdydata);
            C = repmat(normvecs(:,3),1,N_ctrs).*dzrbf(ep,DM_bdydata,dz_bdydata);
            BCM = A + B + C;
            
            % Compute known-terms vector (a.k.a. righthand side vector)
            % This requires the gradient of the unbounded potential at boundary
            gradphi_F_neu = gradphiF_dip(bdydata, srcpnts, dipmom, sig);
            rhs_bdy = -sum(normvecs.*gradphi_F_neu,2);
            
            % Create the full linear system from the blocksand solve it
            % Compose rhs
            rhs = [rhs_int;rhs_bdy];
            % Compose collocation matrix in same order as rhs
            CM = [LCM;BCM];
            
            % Coefficients for evaluation
            %         [coefs,recip_cond] = linsolve(CM,rhs);
            [U,S,V] = svd(CM,0);
            sing_val = diag(S);
            inv_sing_val_TSVD = 1./sing_val;
            indices = find(sing_val < eps);
            inv_sing_val_TSVD(indices) = 0;
            % inv_sing_val_TSVD(floor(0.2*length(inv_sing_val)):end) = 0;
            cbeta = U'*rhs;
            coefs = V*(inv_sing_val_TSVD.*cbeta(1:length(sing_val)));
            
            % Potential at evalpnts in the source free case
            phi0 = EM * coefs;
            % Potential at evalpnts (superposition of effects)
            phi = phi0 + phi_F;
            
            % Compute the difference of the solution with a reference point,
            % arbitrarily chosen as evalpnts(1)
            phi_comp_kansa(:,k) = phi - phi(1);
            elapsed_t_kansa(k) = tt + toc;
            
            % Compute the total errors
            errvec_kansa(k) = errcompute(phi_comp_kansa(:,k),phi_true);
            
            % Store the condition of the system
            condvec_kansa(k) = norm(rhs)/min(sing_val)/...
                norm(cbeta(1:length(sing_val))./sing_val);
            condstr = 'Effective condition number';
            
            Nvec_true_kansa(k) = N_tot;
            
            if iter_out
                fprintf('\terr = %g\n\tcond = %g\n\tN = %d\n',errvec_kansa(k),condvec_kansa(k),Nvec_true_kansa(k));
            end
        end
    end
end

%% Save results in a mat file
save(resultsfilename,'R','sig','evalpnts','dipmom','srcpnts',...
    'phi_comp_*','phi_true');
if BEM
    save(resultsfilename,'*_BEM','-append');
end
if MFS
    save(resultsfilename,'mfs_sphere','mfs_frac','*_MFS','-append');
end
if kansa
    save(resultsfilename,'ep','*_kansa','-append');
end
    
%% Plot simulation results
% Error vs No. of collocation points or elements
h1 = figure;
axes1 = axes('Parent',h1,'FontSize',16,'YScale','log','XLim',[100 6000]);
hold(axes1,'all');
if BEM
    convtest_BEM = semilogy(Nvec_true_BEM, errvec_BEM,...
                            'Parent',axes1,...
                            'MarkerSize',8,...
                            'LineWidth',3);
    set(convtest_BEM,...
        'Marker','hexagram',...
        'Color',[0 0 0],...
        'DisplayName','Symmetric BEM');
end
if kansa
    convtest_kansa = semilogy(Nvec_true_kansa, errvec_kansa,...
                              'Parent',axes1,...
                              'MarkerSize',8,...
                              'LineWidth',3);
    set(convtest_kansa,...
        'Marker','diamond',...
        'Color',[0 0 1],...
        'DisplayName','Kansa (IMQ)');
end
if MFS
    convtest_MFS = semilogy(Nvec_true_MFS,errvec_MFS,...
                            'Parent',axes1,...
                            'MarkerSize',8,...
                            'LineWidth',3);
    set(convtest_MFS(1),...
        'Marker','square',...
        'Color',[1 0 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(1))));
    set(convtest_MFS(2),...
        'Marker','o',...
        'Color',[0 1 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(2))));
    set(convtest_MFS(3),...
        'Marker','*',...
        'Color',[1 0 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(3))));
    set(convtest_MFS(4),...
        'Marker','^',...
        'Color',[0 1 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(4))));
end
title('Current Dipole in Homogeneus Sphere - Convergence test',...
      'FontWeight','bold',...
      'FontSize',20)
xlabel('Total Collocation Points or Elements',...
       'FontWeight','bold',...
       'FontSize',16);
ylabel('Relative Error on Surface Potential',...
       'FontWeight','bold','FontSize',16);
legend1 = legend(axes1,'show');
set(legend1,...
    'Location','SouthWest',...
    'FontWeight','bold',...
    'FontSize',14);
grid on

% CPU time vs No. of collocation points or elements
h2 = figure;
axes1 = axes('Parent',h2,'FontSize',16,'YScale','log','XLim',[100 6000]);
hold(axes1,'all');
if BEM
    time_BEM = semilogy(Nvec_true_BEM, elapsed_t_BEM,...
                        'Parent',axes1,...
                        'MarkerSize',8,...
                        'LineWidth',3);...
    set(time_BEM,...
        'Marker','hexagram',...
        'Color',[0 0 0],...
        'DisplayName','Symmetric BEM');
end
if kansa
    time_kansa = semilogy(Nvec_true_kansa, elapsed_t_kansa,...
                          'Parent',axes1,...
                          'MarkerSize',8,...
                          'LineWidth',3);
    set(time_kansa,...
                          'Marker','diamond',...
                          'Color',[0 0 1],...
                          'DisplayName','Kansa (IMQ)');
end
if MFS
    time_MFS = semilogy(Nvec_true_MFS,elapsed_t_MFS,...
                        'Parent',axes1,...
                        'MarkerSize',8,...
                        'LineWidth',3);
    set(time_MFS(1),...
        'Marker','square',...
        'Color',[1 0 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(1))));
    set(time_MFS(2),...
        'Marker','o',...
        'Color',[0 1 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(2))));
    set(time_MFS(3),...
        'Marker','*',...
        'Color',[1 0 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(3))));
    set(time_MFS(4),...
        'Marker','^',...
        'Color',[0 1 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(4))));
end
title('Current Dipole in Homogeneus Sphere - CPU time',...
      'FontWeight','bold',...
      'FontSize',20)
xlabel('Total Collocation Points or Elements',...
       'FontWeight','bold',...
       'FontSize',16);
ylabel('CPU Time [s]',...
       'FontWeight','bold',...
       'FontSize',16);
legend1 = legend(axes1,'show');
set(legend1,...
    'Location','SouthEast',...
    'FontWeight','bold',...
    'FontSize',14);
grid on

% Effective condition number vs No. of collocation points or elements
h3 = figure;
axes1 = axes('Parent',h3,'FontSize',16,'YScale','log','XLim',[100 6000]);
hold(axes1,'all');
if kansa
    cond_kansa = semilogy(Nvec_true_kansa, condvec_kansa,...
                          'Parent',axes1,...
                          'MarkerSize',8,...
                          'LineWidth',3);
    set(cond_kansa,...
        'Marker','diamond',...
        'Color',[0 0 1],...
        'DisplayName','Kansa (IMQ)');
end
if MFS
    cond_MFS = semilogy(Nvec_true_MFS,condvec_MFS,...
                        'Parent',axes1,...
                        'MarkerSize',8,...
                        'LineWidth',3);
    set(cond_MFS(1),...
        'Marker','square',...
        'Color',[1 0 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(1))));
    set(cond_MFS(2),...
        'Marker','o',...
        'Color',[0 1 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(2))));
    set(cond_MFS(3),...
        'Marker','*',...
        'Color',[1 0 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(3))));
    set(cond_MFS(4),...
        'Marker','^',...
        'Color',[0 1 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(4))));
end
title('Current Dipole in Homogeneus Sphere - System Condition',...
      'FontWeight','bold',...
      'FontSize',20)
xlabel('Total Collocation Points',...
       'FontWeight','bold',...
       'FontSize',16);
ylabel('Effective Condition Number',...
       'FontWeight','bold',...
       'FontSize',16);
legend1 = legend(axes1,'show');
set(legend1,...
    'Location','SouthEast',...
    'FontWeight','bold',...
    'FontSize',14);
grid on

% CPU time vs Error
h4 = figure;
axes1 = axes('Parent',h4,'FontSize',16,'XScale','log','YScale','log');
hold(axes1,'all');
if BEM
    timeerr_BEM = semilogy(errvec_BEM, elapsed_t_BEM,...
                           'Parent',axes1,...
                           'MarkerSize',8,...
                           'LineWidth',3);
    set(timeerr_BEM,...
        'Marker','hexagram',...
        'Color',[0 0 0],...
        'DisplayName','Symmetric BEM');
end
if kansa
    timeerr_kansa = semilogy(errvec_kansa, elapsed_t_kansa,...
                             'Parent',axes1,...
                             'MarkerSize',8,...
                             'LineWidth',3);
    set(timeerr_kansa,...
        'Marker','diamond',...
        'Color',[0 0 1],...
        'DisplayName','Kansa (IMQ)');
end
if MFS
    timeerr_MFS = semilogy(errvec_MFS,elapsed_t_MFS,...
                           'Parent',axes1,...
                           'MarkerSize',8,...
                           'LineWidth',3);
    set(timeerr_MFS(1),...
        'Marker','square',...
        'Color',[1 0 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(1))));
    set(timeerr_MFS(2),...
        'Marker','o',...
        'Color',[0 1 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(2))));
    set(timeerr_MFS(3),...
        'Marker','*',...
        'Color',[1 0 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(3))));
    set(timeerr_MFS(4),...
        'Marker','^',...
        'Color',[0 1 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(4))));
end
title('Current Dipole in Homogeneus Sphere - Cost per Accuracy',...
      'FontWeight','bold',...
      'FontSize',20)
xlabel('Relative Error on Surface Potential',...
       'FontWeight','bold',...
       'FontSize',16);
ylabel('CPU Time [s]',...
       'FontWeight','bold',...
       'FontSize',16);
legend1 = legend(axes1,'show');
set(legend1,...
    'Location','SouthWest',...
    'FontWeight','bold',...
    'FontSize',14);
grid on

% I requested, plot also the analitycal solution, the computed solution and 
% the error distribution for a specific No. of collocation points or 
% elements
if plot_sol
    switch sol_err_style
        case 1
            if BEM
                sol_err_BEM = abs(phi_true - phi_comp_BEM(:,N_plotsol));
                plotstr_BEM = 'BEM Absolute Error [V]';
            end
            if MFS
                sol_err_MFS = abs(phi_true - phi_comp_MFS(:,N_plotsol,mfs_frac_plotsol));
                plotstr_MFS = 'MFS Absolute Error [V]';
            end
            if kansa
                sol_err_kansa = abs(phi_true - phi_comp_kansa(:,N_plotsol));
                plotstr_kansa = 'KM Absolute Error [V]';
            end
        case 2
            if BEM
                sol_err_BEM = log10(abs(phi_true - phi_comp_BEM(:,N_plotsol)));
                plotstr_BEM = 'BEM Log10 of Absolute Error';
            end
            if MFS
                sol_err_MFS = log10(abs(phi_true - phi_comp_MFS(:,N_plotsol,mfs_frac_plotsol)));
                plotstr_MFS = 'MFS Log10 of Absolute Error';
            end
            if kansa
                sol_err_kansa = log10(abs(phi_true - phi_comp_kansa(:,N_plotsol)));
                plotstr_kansa = 'KM Log10 of Absolute Error';
            end
        case 3
            if BEM
                sol_err_BEM = log10(abs(phi_true - phi_comp_BEM(:,N_plotsol))./...
                                    (abs(phi_true)+eps)+eps);
                plotstr_BEM = 'BEM Log10 of Pointwise Relative Error';
            end
            if MFS
                sol_err_MFS = log10(abs(phi_true - phi_comp_MFS(:,N_plotsol,mfs_frac_plotsol))./...
                                    (abs(phi_true)+eps)+eps);
                plotstr_MFS = 'MFS Log10 of Pointwise Relative Error';
            end
            if kansa
                sol_err_kansa = log10(abs(phi_true - phi_comp_kansa(:,N_plotsol))./...
                                      (abs(phi_true)+eps)+eps);
                plotstr_kansa = 'KM Log10 of Pointwise Relative Error';
            end
        otherwise
            error('Unknown 3D plot error style %g',sol_err_style)
    end
    
    h5 = figure; % Analytic
    SurfacePlot_dip(evalpnts, phi_true)
    title('Analytic Potential [V]','FontWeight','bold','FontSize',16)
    if BEM
        h6 = figure;
        subplot(1,2,1)
        SurfacePlot_dip(evalpnts, phi_comp_BEM(:,N_plotsol));
        title('BEM Computed Potential [V]','FontWeight','bold','FontSize',16)
        subplot(1,2,2)
        SurfacePlot_dip(evalpnts, sol_err_BEM);
        title(plotstr_BEM,'FontWeight','bold','FontSize',16)
    end
    if MFS
        h7 = figure;
        subplot(1,2,1)
        SurfacePlot_dip(evalpnts, phi_comp_MFS(:,N_plotsol,mfs_frac_plotsol))
        title('MFS Computed Potential [V]','FontWeight','bold','FontSize',16)
        subplot(1,2,2)
        SurfacePlot_dip(evalpnts, sol_err_MFS);
        title(plotstr_MFS,'FontWeight','bold','FontSize',16)
    end
    if kansa
        h8 = figure;
        subplot(1,2,1)
        SurfacePlot_dip(evalpnts, phi_comp_kansa(:,N_plotsol));
        title('Kansa Computed Potential [V]','FontWeight','bold','FontSize',16)
        subplot(1,2,2)
        SurfacePlot_dip(evalpnts, sol_err_kansa)
        title(plotstr_kansa,'FontWeight','bold','FontSize',16)
    end
end
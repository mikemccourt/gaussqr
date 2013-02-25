% MFS_deriv_test.m
%
% In this test we are interested in studying the evaluation of derivatives
% with the Methods of Fundamental Solutions for a Laplace BVP in 3D.
% We consider a Laplace BVP on a sphere, with Neumann BCs. 
% At the boundary we impose the normal derivatives to be equal to one.
%
% This script allows us to test the convergence rate for increasing N
% (number of collocation points).
%
%  Test parameters:
%     R - Sphere radius <default = 1>
%     mfs_frac - How many centers for MFS, in [0.0,1.0]*N <default = 1.0>
%     mfs_sphere - Fraction beyond R (eg, 1.3R) for centers <default = 1.3>
%
%     Nvec - Row vector of N values <default = 1000:100:1500>
%     N_eval - # points to evaluate error <default = 1001>
%
%  Output parameters:
%     sol_err_style - How do you want the 3D solution error displayed
%                     0 : No error computed, just the solution
%                     1 : Absolute error <default>
%                     2 : Log absolute error
%                     3 : Log pointwise relative error
%     errcolor - Color for error line in log-log plot <default = 'b'>
%     condcolor - Color for condition line in log-log plot <default = 'r'>
%
%  The results of these experiments are stored in
%     err_u_n - Errors computed at Nvec_true
%     cond_CM - Collocation matrix condition numbers at Nvec_true

R = 1;

mfs_frac = 1;
mfs_sphere = 1.2;

N_eval = 1000;
evalpnts = SphereSurfGoldPoints(N_eval,R);

Nvec = 100:400:4100;

sol_err_style = 0;
errcolor = 'b';
condcolor = 'r';

% The normal derivative of the solution is equal to 1 on the boundary
u_n_true = zeros(N_eval,1);
u_b = @(x) sum(x,2);
u_n_true = u_b(evalpnts);

%%%%%%%%%%%%%%%%%%%%%
% Basic setup stuff independent of this problem

% Consider the standard GQR parameters for the errcompute function
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

radbasfun = 'fundamental_3d';
[rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF(radbasfun);

% Loop through the requested N values
err_u_n = [];
cond_CM = [];
k = 1;
for N = Nvec
    fprintf('k=%d\n',k)
    [POINTS,NORMALS] = BallGeometry(R,N,'mfs');
    bdydata = POINTS.bdy11;
    N_ctrs = N*mfs_frac;
    ctrs = SphereSurfGoldPoints(N_ctrs,mfs_sphere*R);
    normvecs = NORMALS.n11;
    
    DM_bdydata = DistanceMatrix(bdydata,ctrs);
    
    % Find all the necessary difference matrices
    dx_bdydata = DifferenceMatrix(bdydata(:,1),ctrs(:,1));
    dy_bdydata = DifferenceMatrix(bdydata(:,2),ctrs(:,2));
    dz_bdydata = DifferenceMatrix(bdydata(:,3),ctrs(:,3));
    
    % Compute collocation matrix for boundary
    A = repmat(normvecs(:,1),1,N_ctrs).*dxrbf(0,DM_bdydata,dx_bdydata);
    B = repmat(normvecs(:,2),1,N_ctrs).*dyrbf(0,DM_bdydata,dy_bdydata);
    C = repmat(normvecs(:,3),1,N_ctrs).*dzrbf(0,DM_bdydata,dz_bdydata);
    CM = A + B + C;
    
    % Compute known-terms vector (a.k.a. righthand side vector)
    % We impose the normal derivative of the function to be 1 on the
    % boundary
    rhs = zeros(N,1);
rhs = u_b(bdydata);
    
    % Solve the system
    [coefs,recip_cond] = linsolve(CM,rhs);
    
    DM_eval = DistanceMatrix(evalpnts, ctrs);
    EM = rbf(0, DM_eval);
    u = EM * coefs;
    
    % Compute the gradient of the solution at evalpnts
    G_u = gradphi0_dip(0, coefs, evalpnts, ctrs, radbasfun);
    
    % Compute the normal derivatives of the solution at evalpnts
    u_n = sum(evalpnts/R.*G_u, 2);
    
    % Compute error on normal derivatives and matrix condition
    err_u_n(k) = errcompute(u_n, u_n_true);
    cond_CM(k) = 1/recip_cond;
    
    fprintf('\terr = %g\n\tcond = %g\n\tN = %d\n',err_u_n(k),cond_CM(k),Nvec(k));
    
    k = k+1;
end

% Plots
clf reset

switch GAUSSQR_PARAMETERS.ERROR_STYLE
    case 1
        errstr = 'Pointwise Rel Err';
    case 2
        errstr = 'Absolute Error';
    case 3
        errstr = 'Relative Error';
    case 4
        errstr = 'RMS Relative Error';
end

[AX,H1,H2] = plotyy(Nvec,err_u_n,Nvec,cond_CM,@loglog);
xlabel('Total collocation points')
set(AX(1),'Xlim',[Nvec(1),Nvec(end)])
set(get(AX(1),'Ylabel'),'String',errstr)
set(AX(1),'Ycolor',errcolor)
set(AX(2),'Xlim',[Nvec(1),Nvec(end)])
set(get(AX(2),'Ylabel'),'String','Matrix Condition')
set(AX(2),'Ycolor',condcolor)
set(H1,'LineWidth',3,'Color',errcolor)
set(H2,'LineWidth',3,'Color',condcolor)
title('MFS derivatives test')

figure
switch sol_err_style
    case 0
        sol_err = u_n;
        plotstr = 'Computed normal derivative';
    case 1
        sol_err = abs(u_n_true - u_n);
        plotstr = 'Absolute Error';
    case 2
        sol_err = log10(abs(u_n_true - u_n));
        plotstr = 'Log10 of Absolute Error';
    case 3
        sol_err = log10(abs(u_n_true - u_n)./(abs(u_n_true)+eps)+eps);
        plotstr = 'Log10 of Pointwise Relative Error';
    otherwise
        error('Unknown 3D plot error style %g',sol_err_style)
end

subplot(1,2,1)
SurfacePlot_dip(evalpnts, u_n_true)
title('True normal derivative','FontWeight','bold','FontSize',12)
subplot(1,2,2)
SurfacePlot_dip(evalpnts, sol_err);
title(plotstr,'FontWeight','bold','FontSize',12)
% The test function
yf = @(x) franke(x(:, 1), x(:, 2));

M = 2;                                  % space dimension 
n = 25;                                 % number of tracks
m = 80;                                 % number of points on tracks
[N, dsites, yy] = TrackData2D(n, m);       % generate N = n*m track data in 2D
neval = 30;                             % parameter for evaluation points
npu = floor((N / 4) ^ (1 / M));          % parameter for PU centres
xx = linspace(0, 1, n);                   % subdomains in one direction
[X, Y] = meshgrid(xx, yy);                % patches centered at tracks
rbf_aniso = @(r) 1 ./ sqrt(1 + r .^ 2);       % define the RBF
rhs = yf(dsites);                        % function values
puctrs = [X(:) Y(:)];                   % define the PU centres
isotropic = false;                       % true/false: require the PU to be isotropic
strat = 'mle';                        % the parametrization strategy ('loocv' / 'mle' / 'no_opt')
tikhonov = 0;                           % a regularization parameter
verbose = true;                         % true/false: print progress to screen
          
global time1 time2 time3
time1 = 0;
time2 = 0;
time3 = 0;
outer_timer = tic;
[epoints, Pf, radii, epsilon] = PU_op_mike(M, dsites, neval, npu, rbf_aniso, yf, rhs, puctrs, isotropic, strat, tikhonov, verbose);
toc(outer_timer)
fprintf('Times in Cost_ep %g %g %g\n', time1, time2, time3)

% What about if we just do interpolation with a C2 Wendland kernel?
% This is just a sanity check for comparison
% If you do not have the GaussQR repo, this will not work, so feel free to
% comment it out.
%
% rbf = @(r) max(1-r,0).^4.*(4*r+1);
% epvec = [10, 5];
% tic
% K = DistanceMatrix(dsites, dsites, epvec, rbf);
% Keval = DistanceMatrix(epoints, dsites, epvec, rbf);
% yeval = Keval * (K \ rhs);
% 
% exact = yf(epoints);
% maxerr = norm(yeval - exact,inf);
% rms_err = norm(yeval - exact)/sqrt(length(exact));
% toc
% fprintf('RMS error:       %e\n', rms_err);
% fprintf('Maximum error:   %e\n', maxerr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Notes from Mike
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The performance of the optimization requires the right choice of
% parameters, but it seems to me that the LOOCV and MLE parametrization
% tools are not performing that well.  This may be because there is
% insufficient "global" behavior when parametrizing locally.  Or maybe
% something else.
%
% You can change the parametrization strategy above, as well as add a
% Tikhonov regularization or choose to run isotropically rather than
% anisotropically.  What I find most confusing (or perhaps interesting) is
% that when I choose the 'no_opt' parametrization strategy, the accuracy of
% the PUM interpolation is very good.  In that setting, the initial guess
% (from the original code) is used on all PU domains with no optimization
% taking place at all.
%
% I am worried that this may have been what was taking place in the
% original experiments in the Dropbox.  When I actually look at the
% optimization taking place in the original code, the LOOCV numbers seem
% kinda bad, especially by comparison to the numbers that we see here.
%
% I think that the problem was the use of the number 25 to denote a bad
% result.  In some of the circumstances, the domain was actually only
% providing values higher than that, leading to the "failure points" being
% the points returned from the optimization.  Try and run the original code
% and look at the LOOCV values, and compare them to the values here to see
% if this actually is consistent with what is computed on your machines.
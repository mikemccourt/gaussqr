%-------------------------------------------------------------------------%
%
% File: PU_for(M,dsites,neval,npu,rbf,wf,f,rhs,r_min,h,P1,ep)
%
% Goal: script that performs partition of unity with variable patches
%       and shape parameters
%
% Inputs:     M:          space dimension
%             dsites:     NXM matrix representing a set of N data sites
%             neval:      number of evaluation points in one direction
%             npu:        number of PU subdomains in one direction
%             rbf:        radial basis function
%             wf:         weight function
%             f:          test function
%             rhs:        function values
%             r_min:      minimum number of points in each PU subdomain
%             h:          upper bound for the radius
%             P1:         parameter for the partition of the PU radius
%             ep:         guess for the shape parameters
%
% Outputs:    epoints:    evaluation points
%             Pf:         interpolant computed at the evaluation points
%
% Calls on:   IntegerBased_MD_Structure, IntegerBased_MD_Neighbourhood,
%             IntegerBased_MD_RangeSearchAniso, MakeSDGrid,
%             IntegerBased_MD_ContainingQuery, DistanceMatrixAniso,
%             PuweightAniso, Cost_ep
%
% Remarks: 1) DistanceMatrixAniso, MakeSDGrid come from the books:
%             [G.E. Fasshauer, Meshfree Approximation Methods with Matlab,
%             World Scientific, Singapore, 2007]
%             [G.E. Fasshauer, M.J. McCourt, Kernel-based Approximation 
%             Methods using Matlab, World Scientific, Singapore, 2015]
%          2) IntegerBased_MD_Structure, IntegerBased_MD_Neighbourhood,
%             IntegerBased_MD_ContainingQuery come from the paper:
%             [R. Cavoretto, A. De Rossi, E. Perracchione, Optimal 
%             selection of local approximants in RBF-PU interpolation,
%             J. Sci. Comput. (2017), in press. 
%             DOI: 10.1007/s10915-017-0418-7]
%
%             fminunc is a Matlab routine used to perform uncostrained 
%             optimization
%
%-------------------------------------------------------------------------%
function [epoints Pf] = PU_for(M,dsites,neval,npu,rbf,...
    f,rhs,r_min,h,P1,ep,puctrs)
% Create neval^M equally spaced evaluation points
epoints = MakeSDGrid(M,neval);
puradius = ones(1,M)*(1./npu);  % Define the initial PU radius
npu_M = size(puctrs,1); neval_M = size(epoints,1); % Initialize;
rbfctrs = dsites; % Define the RBF centres
Pf = zeros(neval_M,1);  % Initialize
% Parameter for integer-based partitioning structure
q = ceil(1./puradius(1,1));
% Build the partitioning structure for data sites and evaluation points
idx_ds = IntegerBased_MD_Structure(dsites,q,puradius(1,1),M);
idx_ep = IntegerBased_MD_Structure(epoints,q,puradius(1,1),M);
% Build the partitioning structure for data sites and evaluation points
radii = zeros(npu_M,M); rn = zeros(npu_M,1); 
epsilon = zeros(npu_M,M); 
for j = 1:npu_M
    puradius = ones(1,M)*(1./npu); old_puradius = puradius; % Initialize
    % Find the block containing the j-th subdomain centre
    index1 = IntegerBased_MD_ContainingQuery(puctrs(j,:),q,...
        puradius(1,1),M);
    % Find data sites located in the j-th subdomain
    [dxx dx] = IntegerBased_MD_Neighbourhood(dsites,idx_ds,index1,q,M,1);
    idx = IntegerBased_MD_RangeSearchAniso(puctrs(j,:),puradius,dxx,dx);
    % Find a range for the PU radius
    t = 1; % Initialize
    while (length(idx) < r_min)
        puradius(:) = puradius(:) + 1/8.*old_puradius(:);
        if max(puradius(:)) > t*max(old_puradius(:))
            t = t + 1;
            [dxx dx] = IntegerBased_MD_Neighbourhood(dsites,idx_ds,...
                index1,q,M,t);
        end
        % Find data sites located in the j-th subdomain
        idx = IntegerBased_MD_RangeSearchAniso(puctrs(j,:),puradius,...
            dxx,dx);
    end
    % Define vectors of radii (semi-axes) to select the optimal ones 
    puradiusvec = zeros(M,P1);
    for d = 1:M
        puradiusvec(d,:) = linspace(puradius(d),h*puradius(d),P1);
    end
    % The LOOCV scheme
    nb = []; ww = 1;
    E = [];
    % NOTE: though all code structure is for MD, the double use of 
    % for-loop makes the code suited to 2D only
    for P = 1:P1
        for Q = 1:P1
            n = 1;
            maxpuradius = max(puradiusvec(1,P),puradiusvec(2,Q));
            while n < 99999
                if maxpuradius > n*max(old_puradius)
                    n = n + 1;
                else
                    t = n;
                    [dxx dx] = IntegerBased_MD_Neighbourhood(dsites,...
                        idx_ds,index1,q,M,t);
                    break
                end
            end
            nb(P,Q) = t;
            % Find data sites located in the j-th subdomain
            idx = IntegerBased_MD_RangeSearchAniso(puctrs(j,:),...
                [puradiusvec(1,P),puradiusvec(2,Q)],dxx,dx);
            % Find optimal ep
            options = optimset('Display','off');
            [minval minerr] = fminunc(@(eparr) Cost_ep(rbf,exp(eparr),...
                idx,dsites,rhs),log(ep),options);            
            epopt = exp(minval);
            minerr = exp(minerr);
            E{ww} = [epopt(1) epopt(2) puradiusvec(1,P) puradiusvec(2,Q)...
                nb(P,Q) minerr];
            ww = ww + 1;
        end
    end
    % Select the radius and the shape parameter
    minE = [];
    for t = 1:ww-1
        minE(t) = [E{1,t}(1,end)];
    end
    [eror ix] = min(minE(:));
    epsilon(j,:) = [E{ix}(1,1) E{ix}(1,2)];
    radii(j,:)= [E{ix}(1,3) E{ix}(1,4)];
    rn(j) = E{ix}(1,5);
    % Construct anisotropic weights
    % Find data sites located in the j-th subdomain
    [dxx dx] = IntegerBased_MD_Neighbourhood(dsites,idx_ds,index1,...
        q,M,rn(j));
    idx = IntegerBased_MD_RangeSearchAniso(puctrs(j,:),radii(j,:),dxx,dx);
    locpts(j).ind = idx;
    [edxx edx] = IntegerBased_MD_Neighbourhood(epoints,idx_ep,index1,...
        q,M,rn(j));
    eidx = IntegerBased_MD_RangeSearchAniso(puctrs(j,:),radii(j,:),...
        edxx,edx);
    elocpts(j).ind = eidx;
end
% Define the pu weight
epu = PuweightAniso(epoints,elocpts,puctrs,radii);
for j = 1:npu_M
    % Compute the distance matrix
    DM_data = DistanceMatrixAniso(dsites(locpts(j).ind,:),...
        rbfctrs(locpts(j).ind,:),epsilon(j,:));
    % Compute the interpolation matrix
    IM = rbf(DM_data);
    if (~isempty(elocpts(j).ind))
        % Compute local evaluation matrix and the local RBF interpolant
        DM_eval = DistanceMatrixAniso(epoints(elocpts(j).ind,:),...
            rbfctrs(locpts(j).ind,:),epsilon(j,:));
        EM = rbf(DM_eval); localfit = EM*(IM\rhs(locpts(j).ind));
        % Accumulate global fit
        Pf(elocpts(j).ind) = Pf(elocpts(j).ind) + localfit.*epu(j).w;
    end
end
% Compute exact solution
exact = f(epoints);
% Compute errors on evaluation grid
maxerr = norm(Pf - exact,inf);
rms_err = norm(Pf - exact)/sqrt(neval_M);
fprintf('RMS error:       %e\n', rms_err);
fprintf('Maximum error:   %e\n', maxerr);
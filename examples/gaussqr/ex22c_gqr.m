% ex22c_gqr
% This is a basic example which demonstrates how RBF-FD works
% In essence, all that is required is a bunch of differentiation matrices
% which get compiled into one big FD operator matrix.
% In the simplest form, some neighboring points are used to approximate the
% derivative at a nearby point.
% Each point gets its own differentiation matrix, which is easy to do in a
% loop.  I'm working on how to do this in a vectorized structure.
%
% NOTE: I may want to introduce a "rescaling" function here to move the
% points to [-1,1]^2 and then account for the scale in the derivative
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = 25;
GAUSSQR_PARAMETERS.ERROR_STYLE = 2; % Need absolute error for 0 solution

% In this example we are testing the quality of the approximation to the
% Laplacian on a tensor product uniform grid
% Easy functions for testing purposes would be harmonic functions, so that
% Lap(f)=0, so we'll test that first.
% Note that the error will only be tested on the interior of the domain,
% since those are the points that have the polynomial-type centered stencil
test_function = @(x) exp(x(:,1)).*sin(x(:,2));
true_solution = @(x) zeros(size(x,1),1);

% Choose the Gaussian shape parameter
% Note that choosing ep small should recover the polynomial result.
ep = 1e-5;

% Choose the number of points in each dimension to use
% Remember that our domain is a tensor product uniform grid [0,1]^2
% We want to consider increasingly dense domains so that h->0
Nvec = 4:15;

% Choose the stencil sizes you want to consider
% Make sure you don't ask for a larger stencil than Nvec can handle
% A is a cell array for the RBF-FD matrices
stencil_size = [5,9,13];
A = cell(length(stencil_size),1);

% Set up error vectors for 5 and 9 point stencils
errvec = cell(length(stencil_size),1);
for k=1:length(errvec)
    errvec{k} = zeros(size(Nvec));
end

h_waitbar = waitbar(0,'Initializing','Visible','on');
    
% Loop through the increasingly dense domain and record the RBF-FD errors
m = 1;
for N=Nvec
    progress = floor(100*((m-1)/range(Nvec)))/100;
    waitbar(progress,h_waitbar,sprintf('%d^2 points in progress',N))
    
    % Set up the interior and whole domain for the problem
    x = pick2Dpoints([0 0],[1 1],N);
    x_int = x(not(any(x==0 | x==1,2)),:);
    
    % Set up the matrices that will be used to compute the Laplacian
    % Note that this matrix will only evaluate the Laplacian on the
    % interior of the domain to match centered polynomial results
    % These are sparse matrices, but we're not working with big problems so
    % I'm not super worried.  This is in here in case we want to worry
    % about preallocation in the future
    for k = 1:length(A)
        A{k} = zeros(size(x_int,1),size(x,1));
    end
    
    % Loop through the interior and find the coefficients for each point
    for x_center = x_int'
        % Sort the points to find the nearest neighbors of xcenter (this_x)
        this_x = x_center';
        [~,sorted_indices] = sort(DistanceMatrix(x,this_x));
        this_x_index_in_x_int = find(DistanceMatrix(x_int,this_x)==0);
        
        % Loop through each of the stencils; compute and store
        for k=1:length(A)
            closest_indices = sorted_indices(1:stencil_size(k))';
            x_closest = x(closest_indices,:);
            
            % Find the Laplacian coefficients for these 5 points
            % This requires some care because of the structured grid
            % This is all that pivoting crap I've been talking about
            GQR = gqr_solveprep(-1,x_closest,ep,1);
            Phi = gqr_phi(GQR,x_closest);
            curr_rank = 1;
            Phi1 = [Phi(:,1),zeros(stencil_size(k),stencil_size(k)-1)];
            phi_ind = 2;
            phi1_ind = 2;
            ind_list = ones(1,stencil_size(k));
            while curr_rank<stencil_size(k)
                Phi1(:,phi1_ind) = Phi(:,phi_ind);
                new_rank = rank(Phi1(:,1:phi1_ind));
                if new_rank>curr_rank
                    ind_list(phi1_ind) = phi_ind;
                    phi1_ind = phi1_ind + 1;
                end
                phi_ind = phi_ind + 1;
                curr_rank = new_rank;
            end
            ind_list = [ind_list,ind_list(end)+1:size(GQR.Marr,2)];
            GQR.Marr = GQR.Marr(:,ind_list);
            GQR.stored_phi1 = Phi1;
            GQR.stored_phi2 = Phi(:,ind_list(stencil_size(k))+1:end);
            invLam1 = diag(1./GQR.eig(GQR.Marr(:,1:stencil_size(k))));
            Lam2 = diag(GQR.eig(GQR.Marr(:,stencil_size(k)+1:end)));
            GQR.Rbar = Lam2*(GQR.stored_phi2'/GQR.stored_phi1')*invLam1;
            psixx = gqr_phi(GQR,this_x,[2 0])*[eye(stencil_size(k));GQR.Rbar];
            psiyy = gqr_phi(GQR,this_x,[0 2])*[eye(stencil_size(k));GQR.Rbar];
            Psi = GQR.stored_phi1 + GQR.stored_phi2*GQR.Rbar;
            Lap_Mat_this_x = (psixx + psiyy)/Psi;
            
            % Plug these values into the matrix
            % We use the fact that the first value in closest_indices is
            % the index of this_x, since it is the closest to itself
            A{k}(this_x_index_in_x_int,closest_indices) = Lap_Mat_this_x;
        end
    end
    
    % Store the data and the true solution
    y = test_function(x);
    yL = true_solution(x_int);
    
    % Compute each approximate Laplacian and the error
    for k=1:length(A)
        yL_appx = A{k}*y;
        errvec{k}(m) = errcompute(yL_appx,yL);
    end
    m = m + 1;pause
end

waitbar(1,h_waitbar,'Plotting')

% Plot the output of these RBF-FD tests
h = figure;
ph = zeros(size(errvec));
legend_strings = cell(length(errvec),1);
hold all
for k=1:length(errvec)
    ph(k) = loglog(Nvec.^2,errvec{k},'linewidth',3);
    pfit = polyfit(2*log(Nvec),log(errvec{k}),1);
    legend_strings{k} = sprintf('width %d, slope %3.2f',stencil_size(k),pfit(1));
    set(get(ph(k),'parent'),'xscale','log')
    set(get(ph(k),'parent'),'yscale','log')
end
hold off
title('Laplacian of a harmonic function')
xlabel('Number of points in [0,1]^2')
ylabel('Absolute Error')
legend(ph,legend_strings,'location','northwest');

close(h_waitbar)
% ex5e_gqr
% This problem compares GaussQR collocation with a very refined finite
% difference solution
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .6;

clear functions % reset persistent variables

% Discretization choices
N = 25;
dt = .1;
T = 10;

% Solution storage choices
t_store = 0:5:T;
sol_store = cell(size(t_store));
count_store = 1;

% Solution parameters
h = 1e-7; % Jacobian-vector Finite difference parameter
gmres_restart = 30;
Nit_max = 20; % Max Newton iterations
du_tol = 1e-5; % Newton step length for convergence

% Set up spatial, temporal discretization
time_steps = dt:dt:T;
x = pick2Dpoints([-1,-1],[1,1],N);

% Line search values
tls = [0 .25 .75 1];
Fls = zeros(1,4);

% User preferences
plot_initial_conditions = 0;
plot_solutions = 0;
keep_me_notified = 1; % 0 - no, 1 - Newton level, 2 - gmres level
store_solutions = 1;

% 0 - only compute the preconditioner once per time step
% 1 - direct solve at each Newton step
recompute_preconditioner = 1; 
% 0 - use the full Jacobian as the preconditioner
% 1 - use the ILUdt as the preconditioner
use_ilu_factorization = 0; % iterative solve using ILU

u = ex5e_gqr_res(x); % Evaluate initial condition

if plot_initial_conditions
    X = reshape(x(:,1),N,N);
    Y = reshape(x(:,2),N,N);
    A = reshape(u(1:N*N),N,N);
    H = reshape(u(N*N+1:end),N,N);
    subplot(1,2,1)
    surf(X,Y,A),title('Activator'),xlabel('x'),ylabel('y')
    subplot(1,2,2)
    surf(X,Y,H),title('Inhibitor'),xlabel('x'),ylabel('y')
    pause
end

if store_solutions
    sol_store{count_store} = u;
    count_store = count_store + 1;
end

% Store computational time
computing_time = zeros(size(time_steps));

% Solve each time step using the Newton method
%   J*du = -Fu
%   u_new = u_old + du
uold = u;
step = 1;
for t=time_steps
    normu = norm(uold);
    tic % Start the timer
    for k=1:Nit_max
        % Create the preconditioner
        if recompute_preconditioner | k==1
            if keep_me_notified==2
                if k==1
                    fprintf('Computing preconditioner for first Newton step\n')
                else
                    fprintf('Recomputing preconditioner, as requested\n')
                end
            end
            J = ex5e_gqr_res(u,x,uold,dt,1);
            if use_ilu_factorization & ~recompute_preconditioner
                [L,U] = luinc(J,1e-6);
            end
        end

        % Evaluate the RHS
        Fu = -ex5e_gqr_res(u,x,uold,dt);

        % Solve the system
        if recompute_preconditioner
            du = J\Fu;
            flag = 0;
            iter = 1;
        elseif use_ilu_factorization
            [du,flag,relres,iter] = gmres(@(b)ex5e_gqr_res(u,x,uold,dt,h,b),Fu,gmres_restart,1e-15,[],L,U);
        else
            [du,flag,relres,iter] = gmres(@(b)ex5e_gqr_res(u,x,uold,dt,h,b),Fu,gmres_restart,[],[],J);
        end
        
        % Check outcome of GMRES, and compute solution norm
        if flag>0
            error('gmres did not converge, flag=%d, iter=%d',flag,iter(1))
        elseif keep_me_notified==2
            fprintf('solver converged in %d outer iterations\n',iter(1))
        end
        normdu = norm(du);
        
        % Perform the cubic line search
        Fls(1) = norm(Fu);
        Fls(2) = norm(ex5e_gqr_res(u+.25*du,x,uold,dt));
        Fls(3) = norm(ex5e_gqr_res(u+.75*du,x,uold,dt));
        Fls(4) = norm(ex5e_gqr_res(u+du,x,uold,dt));
        ls_val = fminbnd(@(t)polyval(polyfit(tls,Fls,3),t),0,1);
        
        % Update the solution to the latest time step
        u = u + ls_val*du;
        
        % Check for convergence or divergence
        if normdu/normu<du_tol
            if keep_me_notified
                fprintf('Newton convergence: time %g, %d steps, %g tol\n',t,k,normdu)
            end
            break
        elseif k==Nit_max
            warning('No Newton convergence in %d steps',Nit_max)
        end
    end
    uold = u;
    computing_time(step) = toc;

    if store_solutions
        if t==t_store(count_store)
            sol_store{count_store} = u;
            count_store = count_store + 1;
        end
    end

    step = step + 1;
end

if plot_solutions
    X = reshape(x(:,1),N,N);
    Y = reshape(x(:,2),N,N);
    A = reshape(u(1:N*N),N,N);
    H = reshape(u(1+N*N:end),N,N);
    subplot(1,2,1)
    surf(X,Y,A),title('Activator'),xlabel('x'),ylabel('y')
    subplot(1,2,2)
    surf(X,Y,H),title('Inhibitor'),xlabel('x'),ylabel('y')
end

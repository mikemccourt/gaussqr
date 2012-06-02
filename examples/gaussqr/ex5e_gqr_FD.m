function ex5e_gqr_FD(N,dt,t_store,fileName)
% function ex5e_gqr_FD(N,dt,t_store,fileName)
% This computes the finite difference solution to the Geirer-Meinhardt
% problem so that we can make a comparison to the GaussQR solution.  It
% will store the values sol_store, t_store, N, x in the file gmFDsols.mat
% Inputs : N - number of FD points in each dimension (N^2 total points)
%          dt - time step size
%          t_store - the times at which you require the solution
%          fileName - name of the file to store to
% Output : Stored in fileName - sol_store
%                               t_store
%                               N
%                               x
%
% The temporal discretization is Crank-Nicolson

% In case an unsorted vector
t_store = sort(t_store);

% Solution storage setup
sol_store = cell(size(t_store));
count_store = 1;

% Solution parameters
h = 1e-7; % Jacobian-vector Finite difference parameter
gmres_restart = 30;
Nit_max = 20; % Max Newton iterations
du_tol = 1e-7; % Newton step length for convergence

% Set up spatial, temporal discretization
time_steps = dt:dt:t_store(end);
x = pick2Dpoints([-1,-1],[1,1],N);

% Line search values
tls = [0 .25 .75 1];
Fls = zeros(1,4);

% -1 - print nothing to the screen
% 0 - only print when a solution is stored
% 1 - print outputs at the Newton solve level
% 2 - print outputs at the linear solve level
keep_me_notified = 0;

% 0 - only compute the preconditioner once per time step
% 1 - direct solve at each Newton step
recompute_preconditioner = 1; 
% 0 - use the full Jacobian as the preconditioner
% 1 - use the ILUdt as the preconditioner
use_ilu_factorization = 0; % iterative solve using ILU

u = ex5e_gqr_FD_res(x); % Evaluate initial condition

if t_store(count_store)==0
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
            J = ex5e_gqr_FD_res(u,x,uold,dt,1);
            if use_ilu_factorization & ~recompute_preconditioner
                [L,U] = luinc(J,1e-6);
            end
        end

        % Evaluate the RHS
        Fu = -ex5e_gqr_FD_res(u,x,uold,dt);

        % Solve the system
        if recompute_preconditioner
            du = J\Fu;
            flag = 0;
            iter = 1;
        elseif use_ilu_factorization
            [du,flag,relres,iter] = gmres(@(b)ex5e_gqr_FD_res(u,x,uold,dt,h,b),Fu,gmres_restart,1e-15,[],L,U);
        else
            [du,flag,relres,iter] = gmres(@(b)ex5e_gqr_FD_res(u,x,uold,dt,h,b),Fu,gmres_restart,[],[],J);
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
        Fls(2) = norm(ex5e_gqr_FD_res(u+.25*du,x,uold,dt));
        Fls(3) = norm(ex5e_gqr_FD_res(u+.75*du,x,uold,dt));
        Fls(4) = norm(ex5e_gqr_FD_res(u+du,x,uold,dt));
        ls_val = fminbnd(@(t)polyval(polyfit(tls,Fls,3),t),0,1);
        
        % Update the solution to the latest time step
        u = u + ls_val*du;
        
        % Check for convergence or divergence
        if normdu/normu<du_tol
            if keep_me_notified
                fprintf('\t Newton convergence: time %g, %d steps, %g tol\n',t,k,normdu)
            end
            break
        elseif k==Nit_max
            warning('No Newton convergence in %d steps',Nit_max)
        end
    end
    uold = u;
    computing_time(step) = toc;

    if abs(t-t_store(count_store))<1e-7
        sol_store{count_store} = u;
        count_store = count_store + 1;
        if keep_me_notified>=0
            fprintf('\t\t Finite difference solution stored at t=%g\n',t)
        end
    end

    step = step + 1;
end

save(fileName,'sol_store','t_store','N','x');
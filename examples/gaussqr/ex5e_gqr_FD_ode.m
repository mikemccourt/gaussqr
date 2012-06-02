function [computing_time,time_steps] = ex5e_gqr_FD_ode(N,dt,t_store,fileName)
% function [computing_time,time_steps] = ex5e_gqr_FD_ode(N,dt,t_store,fileName)
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
% It uses a built-in Matlab ODE solver

% In case an unsorted vector
t_store = sort(t_store);
T = t_store(end);

% Solution storage setup
sol_store = cell(size(t_store));
count_store = 1;

% Set up spatial discretization
x = pick2Dpoints([-1,-1],[1,1],N);

% Set up initial conditions and mass matrix
[u0,M] = ex5e_gqr_FD_ode_rhs(x);

% Set up the ODE solver options
options = odeset('InitialStep',dt ...
                 ,'Jacobian',@(t,u)ex5e_gqr_FD_ode_rhs(x,u) ...
                 ,'Mass',M ...
                 ,'MStateDependence','none' ...
                 ,'MassSingular','yes' ...
                 ,'OutputFcn',@ex5e_gqr_output_function ...
                 );

% Solve the system
sol = ode15s(@(t,u)ex5e_gqr_FD_ode_rhs(x,u,t),[0,T],u0,options);

% Evaluate the solution at the requested points
solmat = deval(sol,t_store);
for k=1:length(t_store)
    sol_store{k} = solmat(:,k);
end

% Save the data to the requested file
save(fileName,'sol_store','t_store','N','x');

% Return the required outputs
time_steps = sol.x;
computing_time = ex5e_gqr_output_function([],[],'time');

end

function status = ex5e_gqr_output_function(t,u,flag)
persistent comp_time

status = 0;
if length(flag)==0
    fprintf('Solution found at t=%g\n',t)
    comp_time = [comp_time,toc];
    tic
else
    switch flag
        case 'init'
            fprintf('Starting time stepping\n')
            tic
        case 'done'
            fprintf('Completing time stepping\n')
            comp_time = [comp_time,toc];
        case 'time'
            status = comp_time;
            comp_time = [];
        otherwise
            fprintf('Unknown flag=%s\n',flag)
            status = 1
    end
end
end

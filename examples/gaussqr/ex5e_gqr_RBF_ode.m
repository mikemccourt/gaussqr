function [computing_time,time_steps] = ex5e_gqr_RBF_ode(N,t_store,fileName)
% function [computing_time,time_steps] = ex5e_gqr_RBF_ode(N,t_store,fileName)
% This computes the finite difference solution to the Geirer-Meinhardt
% problem so that we can make a comparison to the GaussQR solution.  It
% will store the values sol_store, t_store, N, x in the file gmFDsols.mat
% Inputs : N - number of FD points in each dimension (N^2 total points)
%          t_store - the times at which you require the solution
%          fileName - name of the file to store to
% Output : Stored in fileName - sol_store
%                               t_store
%                               N
%                               x
%
% It uses a built-in Matlab ODE solver

% In case an unsorted vector is passed
t_store = sort(t_store);
T = t_store(end);

% Solution storage setup
sol_store = cell(size(t_store));
count_store = 1;

% Set up spatial discretization
x = pick2Dpoints([-1,-1],[1,1],N);

% Set up initial conditions and mass matrix
[u0,M] = ex5e_gqr_RBF_ode_rhs(x);

xx = pick2Dpoints([-1,-1],[1,1],3*N);
uu = ex5e_gqr_RBF_ode_rhs(xx);

% Test ability for GaussQR to reproduce initial conditions
alpha = 1;
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = ceil(.1*N*N);
epvec = logspace(-1,1,25);
errvec = zeros(size(epvec));
y = u0(N*N+1:end);
yy = uu(9*N*N+1:end);
k = 1;
for ep=epvec
    GQR = gqr_rsolve(x,y,ep,alpha);
    yp = gqr_eval(GQR,xx);
    errvec(k) = errcompute(yp,yy);
    k = k+1;
end
loglog(epvec,errvec)
pause

% X = reshape(x(:,1),N,N);
% Y = reshape(x(:,2),N,N);
% A = reshape(u0(1:N*N),N,N);
% H = reshape(u0(1+N*N:end),N,N);
% subplot(1,2,1)
% surf(X,Y,A)
% subplot(1,2,2)
% surf(X,Y,H)
% pause

% Set up the ODE solver options
options = odeset('Mass',M ...
                 ,'MStateDependence','none' ...
                 ,'MassSingular','yes' ...
                 ,'OutputFcn',@ex5e_gqr_output_function ...
                 );

% % Solve the system
% sol = ode15s(@(t,u)ex5e_gqr_RBF_ode_rhs(x,u,t),[0,T],u0,options);

uold = u0;
sol_store{1} = u0;

u = uold;
h = 1e-7; % Finite diff parameter
for t_step=2:length(t_store);
    t = t_store(t_step);
    dt = t-t_store(t_step-1);
    for k=1:20
        rhs = ex5e_gqr_RBF_ode_rhs(x,u,t);
        Fu = M*(u-uold)/dt - rhs;
        delta_u = gmres(@(b)M*b/dt - (ex5e_gqr_RBF_ode_rhs(x,u+h*b,t) - rhs)/h,-Fu);
        u = u + .4*delta_u;
        fprintf('\tnorm(Fu)=%g, norm(delta_u)=%g\n',norm(Fu),norm(delta_u));
    end
    sol_store{t_step} = u;
    uold = u;
    
    X = reshape(x(:,1),N,N);
    Y = reshape(x(:,2),N,N);
    A = reshape(u(1:N*N),N,N);
    H = reshape(u(1+N*N:end),N,N);
    subplot(1,2,1)
    surf(X,Y,A),title(sprintf('Activator, t=%g',t));
    subplot(1,2,2)
    surf(X,Y,H),title(sprintf('Inhibitor, t=%g',t));
    pause
end

% % Evaluate the solution at the requested points
% solmat = deval(sol,t_store);
% for k=1:length(t_store)
%     sol_store{k} = solmat(:,k);
% end

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
    fprintf('RBF solution found at t=%g\n',t)
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

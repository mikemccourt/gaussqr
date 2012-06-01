% ex5e_gqr
% This problem compares GaussQR collocation with a very refined finite
% difference solution
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .6;

clear functions % reset persistent variables

N = 5;
dt = .001;
T = .001;
h = 1e-7;
gmres_restart = 30;
Nit_max = 20;
du_tol = 1e-5;

x = pick2Dpoints([-1,-1],[1,1],N);

% Line search values
tls = [0 .25 .75 1];
Fls = zeros(1,4);

plot_initial_conditions = 0;

uold = ex5e_gqr_res(x); % Evaluate initial condition

if plot_initial_conditions
    N2 = size(x,1);
    aold = uold(1:N2);
    hold = uold(N2+1:end);
    X = reshape(x(:,1),N,N);
    Y = reshape(x(:,2),N,N);
    A = reshape(aold,N,N);
    H = reshape(hold,N,N);
    subplot(1,2,1)
    surf(X,Y,A),title('Activator')
    subplot(1,2,2)
    surf(X,Y,H),title('Inhibitor')
end

% Solve each time step using the Newton method
%   J*du = -Fu
u = uold;
for t=dt:T
    for k=1:Nit_max
        % Create the preconditioner
        J = ex5e_gqr_res(u,x,uold,dt,1);

        % Evaluate the RHS
        Fu = -ex5e_gqr_res(u,x,uold,dt);

        % Solve the system
        du = gmres(@(b)ex5e_gqr_res(u,x,uold,dt,h,b),Fu,gmres_restart,[],[],J);
        
        % Perform the line search
        Fls(1) = norm(Fu);
        Fls(2) = norm(ex5e_gqr_res(u+.25*du,x,uold,dt));
        Fls(3) = norm(ex5e_gqr_res(u+.75*du,x,uold,dt));
        Fls(4) = norm(ex5e_gqr_res(u+du,x,uold,dt));
        ls_val = fminbnd(@(t)polyval(polyfit(tls,Fls,3),t),0,1);
        
        % Update the solution to the latest time step
        u = u + ls_val*du;
        
        if norm(du)<du_tol
            break
        elseif k==Nit_max
            warning('No Newton convergence in %d steps',Nit_max)
        end
    end
    uold = u;
end
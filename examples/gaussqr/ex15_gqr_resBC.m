function fu = ex15_gqr_resBC(coef,GQR,x,uold,dt,BC,t)
% function fu = ex15_gqr_resBC(coef,GQR,x,uold,dt,BC,t)
% This function takes in a set of coefficients c and returns the residual
% evaluated using those coefficients in the eigenfunction series
%   coef - summation coefficients to use in GQR
%   GQR  - GaussQR evaluation object
%   x    - points at which to evaluate residual
%   uold - previous solution
%   dt   - time step
%   BC   - 2-vector of boundary conditions
%   t    - current time to evaluate at

% Need the true solution for boundary conditions
    usol = @(x,t) exp(-t)*(1-x.^2);
    usolx = @(x,t) exp(-t)*(-2*x);
    
% Apply the test coefficients
    GQR.coef = coef;
    
% Not sure that I need this, but I'm keeping it for now
    BCenforce = 1/dt;
%     BCenforce = 1;
    
% This evaluates the residual on the interior
% The residual is F(u) = u_t - (k(u_x)u_x)_x - f
    u = gqr_eval(GQR,x);
    u_x = gqr_eval(GQR,x,1);
    u_xx = gqr_eval(GQR,x,2);
    u_t = (u-uold)/dt;
    ku_x = kfunc(u_x);
    kpu_x = kfunc(u_x,1);
    f_source = ffunc(x,t);
    fu = u_t - u_xx.*(kpu_x.*u_x+ku_x) - f_source;
    
% Apply the boundary conditions
    switch BC(1)
        case 0 % Dirichlet
            fu(1) = gqr_eval(GQR,x(1))-usol(x(1),t);
        case 1 % Neumann
            fu(1) = gqr_eval(GQR,x(1),1)-usolx(x(1),t);
    end
    
    switch BC(end)
        case 0 % Dirichlet
            fu(end) = gqr_eval(GQR,x(end))-usol(x(end),t);
        case 1 % Neumann
            fu(end) = gqr_eval(GQR,x(end),1)-usolx(x(end),t);
    end

% Used to scale the boundary conditions up to the interior
    fu([1,end]) = fu([1,end])*BCenforce;
end

% This is the source term
function f = ffunc(x,t)
    u_x = exp(-t)*(-2*x);
    u_xx = exp(-t)*(-2);
    u_t = -exp(-t)*(1-x.^2);
    ku_x = kfunc(u_x);
    kpu_x = kfunc(u_x,1);
    f = u_t-u_xx.*(kpu_x.*u_x+ku_x);
end

function k = kfunc(u_x,deriv)
% function k = kfunc(u_x,deriv)
% This function evaluated the diffusivity given the gradient of the
% temperature within the reactor
% Note right now that passing deriv=1 will produce the derivative with
% respect to u_x, not the derivative wrt to GQR.coef.  I'm not totally sure
% yet what you need to do to get that Jacobian
%
% Eventually I may need to more carefully evaluate log(cosh) because the
% actual value may be fine, but evaluating it directly may be a problem
%
% Setting kk = 0 will convert the problem back to a linear problem
    kk = 0;
    z = 2;
    C = 1;
    k0 = 1;
    B = .5*kk/z*log(1+cosh(2*z*C))-kk*C+(kk-2)*.5*log(2)/z;
    if nargin==1
        k = .5*kk/z*log(cosh(2*z*u_x)+cosh(2*z*C))-kk*C+(kk-2)*.5*log(2)/z-B+k0;
    else
        k = kk*sinh(2*z*u_x)./(cosh(2*z*u_x)+cosh(2*z*C));
    end
end
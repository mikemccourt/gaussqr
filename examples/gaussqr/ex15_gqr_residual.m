function fu = ex15_gqr_residual(GQR,x,uold,dt,b,h)
% function fu = ex15_gqr_residual(GQR,x,uold,dt)
% This function takes in a set of coefficients c and returns the residual
% evaluated using those coefficients in the eigenfunction series
%   GQR  - GQR.coef has coefficients for residual evaluation
%   x    - points at which to evaluate residual
%   dt   - time step
%   uold - previous solution
%
% function fu = ex15_gqr_residual(GQR,x,uold,dt,b,h)
% When h is passed, it is used to approximate Jacobian-vector products
% involving F(u), as evaluated with the function above
%   b - vector (or set of vectors) to multiply to, J(F)(u)b
%   h - finite differencing parameter for matrix-free Jacobian evaluation
%
% Note that this only works with the PDE, and to enforce boundary
% conditions you'll need to do some post processing
    if nargin==4
        u = gqr_eval(GQR,x);
        u_x = gqr_eval(GQR,x,1);
        u_xx = gqr_eval(GQR,x,2);
        u_t = (u-uold)/dt;
        ku_x = kfunc(u_x);
        kpu_x = kfunc(u_x,1);
        fu = u_t-u_xx.*(kpu_x.*u_x+ku_x);
    elseif nargin==6
        m = size(b,2);
        u = ex15_gqr_residual(GQR,x,uold,dt);
        c = GQR.coef;
        
        fu = zeros(size(b));
        for k=1:m
            GQR.coef = c+h*b(:,k);
            uhb = ex15_gqr_residual(GQR,x,uold,dt);
            fu(:,k) = (uhb-u)/h;
        end
        
        GQR.coef = c; % Don't know if I need to reset it
    else
        error('Wrong number of paramaeters')
    end
end

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
    kk = 1;
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
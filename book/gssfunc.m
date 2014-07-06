function [fx,J] = gssfunc(v,gval)
% function [fx,J] = gssfunc(v,gval)
% Inputs: v - vector to evaluate F(v) at
%         gval - (optional) inner product BC magnifier <default=1>
% Outputs: fx - F(v) value for the passed v
%               NOTE: negative v(5) or v(6) values yield NaN
%          J - J(F)(v) Jacobian evaluated at v
% NOTE: If this function is called with no inputs, it will return a
%       optimset options object full of default values for the fsolve
%       routine and a default initial guess:
%            [options,x0] = gssfunc;
%            gval = 5;
%            [xsol,...] = fsolve(@(v)gssfunc(v,gval),x0,options);
%
% This function is the F(v) function which should allow us to solve for the
% coefficients which define the inner product for the Sobolev space
% associated with the operators
%       P = d/dx
%       B = I_[0,1]
% such that the inner product is
%   (f,g)_{A,C,P,B} = (f,g)_P + gval*(f,g)_B
%
% The reproducing kernel for this inner product is
%   K(x,z) = G(x,z) + R(x,z)
% where
%   (P*P)G = delta
%       BG = 0
% defines G(x,z) = min(x,z) - x*z and
%   (P*P)B = 0
% defines R(x,z) = v(5)*phi_1(x)phi_1(z) + v(6)*phi_2(x)phi_2(z).  The phi
% basis functions are defined through the psi basis functions
%   psi_1(x) = x   &   psi_2(x) = 1-x
% such that
%   phi_1 = v(1)psi_1 + v(2)psi_2   &   phi_2 = v(3)psi_1 + v(4)psi_2
% The psi are designed this way to be orthonormal in the B inner product
%
% The equations we want to solve are the same B-ON equations:
%    (phi_1,phi_2)_B = 0
%    (phi_1,phi_1)_B = 1
%    (phi_2,phi_2)_B = 1
% as well as the equations created by making the inner product from Qi's
% paper match the inner product described above.  That matching yields 4
% equations, but 2 are redundant, so there are 6 equations in 6 unknowns.
% The terms that are matched are f(0)g(0), f(1)g(1) and f(1)g(0).
%
% To determine these coefficients when gval = 5, you run
%   gval = 5;
%   x0 = zeros(6,1);
%   gss_coef = fsolve(@(x)gssfunc(x,gval),x0);

% If you are just calling this to get the default optimization routine
% values by calling this with no arguments passed
if nargin==0
    fx=optimset('LargeScale','off',...
                'Display','iter',...
                'FunValCheck','on',...
                'DerivativeCheck','on',...
                'Diagnostics','on',...
                'Jacobian','on',...
                'NonlEqnAlgorithm','lm',...
                'LineSearchType','cubicpoly');
    J = ones(6,1);
    return
end

% This is the default case
if not(exist('gval','var'))
    gval = 1;
end

% Handle the degenerate cases of 0 values for a_1 & a_2 (v(5) & v(6))
a = v(1);
b = v(2);
c = v(3);
d = v(4);
if v(5)<0
    e = NaN;
elseif v(5)==0
    e = Inf;
else
    e = v(5);
end
if v(6)<0
    f = NaN;
elseif v(6)==0
    f = Inf;
else
    f = v(6);
end

% We're really solving F(v)=y, so we subtract the constant terms in a
% separate vector to help me think about things
fx = [a*c+b*d;
      a^2+b^2;
      c^2+d^2;
      b^2/e+d^2/f-(1-2*a*b)*b^2-(1-2*c*d)*d^2+2*(a*d+b*c)*b*d;
      a^2/e+c^2/f-(1-2*a*b)*a^2-(1-2*c*d)*c^2+2*(a*d+b*c)*a*c;
      a*b/e+c*d/f-(1-2*a*b)*a*b-(1-2*c*d)*c*d+(a*d+b*c)^2] ...
      - [0;1;1;gval;gval;0];

% Evaluate the Jacobian as needed
% I'm not sure how I feel about doing this given the discontinuous
% derivative in the a_1 & a_2 coefficients (not being present when those
% values are 0) but we'll see if this works.
if nargout>1
    J = [c,d,a,b,0,0;
         2*a,2*b,0,0,0,0;
         0,0,2*c,2*d,0,0;
         2*b^3+2*b*d^2,2*b/e+6*a*b^2-2*b+4*a*b*d,2*d^3+2*b^2*d,2*d/f+6*c*d^2-2*d+4*a*b*d,-b^2/e^2,-d^2/f^2;
         2*a/e+6*a^2*b-2*a+4*a*c*d,2*a^3+2*a^2*c,2*c/f+6*c^2*d-2*c+4*a*b*c,2*c^3+2*a^2*c,-a^2/e^2,-c^2/f^2;
         b/3+4*a*b^2-b+2*(a*d+b*c)*d,a/3+4*a^2*b-a+2*(a*d+b*c)*c,d/f+4*c*d^2-d+2*(a*d+b*c)*b,c/f+4*c^2*d-c+2*(a*d+b*c)*a,-a*b/e^2,-c*d/f^2];
end
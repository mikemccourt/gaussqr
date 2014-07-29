function [fx,J] = gssfunc(v,ipcoef)
% function [fx,J] = gssfunc(v,ipcoef)
% Inputs: v - vector to evaluate F(v) at
%         ipcoef - inner product BC coefficients
% Outputs: fx - F(v) value for the passed v
%               NOTE: negative v(3) or v(4) values yield NaN
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
%   (f,g)_{A,C,P,B} = (f,g)_P + ipcoef(1)*f0*g0
%                             + ipcoef(2)*(f0*g1+f1*g0)
%                             + ipcoef(3)*f1*g1
% We have made the assumption that the definition of the boundary inner
% product is:
%   (f,g)_partialOmega = f0*g0 + f1*g1
% although we may consider a more general definition evenutally.
%
% The reproducing kernel for this inner product is
%   K(x,z) = G(x,z) + R(x,z)
% where
%   (P*P)G = delta
%       BG = 0
% defines G(x,z) = min(x,z) - x*z and
%   (P*P)B = 0
% defines R(x,z) = v(3)*phi_1(x)phi_1(z) + v(4)*phi_2(x)phi_2(z).  The phi
% basis functions are defined through the psi basis functions
%   psi_1(x) = x   &   psi_2(x) = 1-x
% such that
%   phi_1 = v(1)psi_1 + v(2)psi_2   &   phi_2 = -/+v(2)psi_1 +/- v(1)psi_2
% The psi are designed this way to be orthonormal in the B inner product.
% The +/- ambiguity in phi_2 goes away in K.
%
% The equations we want to solve are the same B-ON equations:
%    (phi_1,phi_2)_B = 0
%    (phi_1,phi_1)_B = 1
%    (phi_2,phi_2)_B = 1
% Some moving content around and cleaning things up allows to enforce only
% (phi_1,phi_1)_B = 1.  After the user chooses alpha, beta, gamma from
% above we define the full set of 4 equations as
%                    (1-v(1)^2)/v(3) + v(1)^2/v(4) - 1 = ipcoef(1)
%                        v(1)*v(2)*(1/v(3)-1/v(4)) + 1 = ipcoef(2)
%                    v(1)^2/v(3) + (1-v(1)^2)/v(4) - 1 = ipcoef(3)
%                                      v(1)^2 + v(2)^2 = 1
%
% Example:
% To determine these coefficients when the inner product you want is
%      (f,g)_{A,C,P,B} = (f,g)_P + 4*f0*g0 + 2*(f0*g1+f1*g0) + 2*f1*g1
%   ipcoef = [3 -1 1];
%   [gss_opts,x0] = gssfunc; % Get the optimization options
%   gss_coef = fsolve(@(x)gssfunc(x,ipcoef),x0,gss_opts);

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
    J = ones(4,1);
    return
elseif nargin~=2
    error('Incorrect number of arguments passed, nargin=%d',nargin)
end

% Here we store the coefficients in simple form
% The solution requires a1,a2>0, so we allow the solver to explore negative
% a1, a2 but always use the absolute value in the function evaluation and
% in defining the kernel
c1 = v(1);
c2 = v(2);
a1 = abs(v(3));a2 = abs(v(4));
s1 = sign(v(3));s2 = sign(v(4));
alpha = ipcoef(1);
beta = ipcoef(2);
gamma = ipcoef(3);
% Eventually, I will need to
% handle the degenerate cases of 0 values for a_1 & a_2 (v(3) & v(4))
% I think this will have to be more complicated now since there are
% different equations when a_1 or a_2 is zero

% We're really solving F(v)=y, so we subtract the constant terms in a
% separate vector to help me think about things
fx = [(1-c1^2)/a1 + c1^2/a2 - 1
      c1*c2*(1/a1 - 1/a2) + 1
      c1^2/a1 + (1-c1^2)/a2 - 1
      c1^2 + c2^2] ...
      - [alpha;beta;gamma;1];

% Evaluate the Jacobian as needed
% I'm not sure how I feel about doing this given the discontinuous
% derivative in the a_1 & a_2 coefficients (not being present when those
% values are 0) but I guess this is okay for now.
if nargout>1
    J = [-2*c1/a1+2*c1/a2, 0, (c1^2-1)/a1^2, -c1^2/a2^2
         c2*(1/a1-1/a2), c1*(1/a1-1/a2), -c1*c2/a1^2, c1*c2/a2^2
         2*c1/a1 - 2*c1/a2, 0, -c1^2/a1^2, (c1^2-1)/a2^2
         2*c1, 2*c2, 0, 0];
    J = J*diag([1,1,s1,s2]); % To account for the absolute value in the a1, a2 def above
end
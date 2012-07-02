function [ g ] = genBoundarySatisfyingFunction( f, beta, lbound, ubound, lbval, ubval )
%GENBOUNDARYSATISFYINGFUNCTION
% Takes a symbolic f(x) and generates the function g(x) = f(x)+(terms
%that cause f and all even derivatives of f up to the (2*beta)th derivative
%to equal 0 at the at the boundary of [lbound,ubound])
%   Returns a function g based on f that satisfies the boundary 
%   conditions
%      f(lbound) = f''(lbound) = f''''(lbound) = f^(2*beta)(lbound) = lcond
%      f(ubound) = f''(ubound) = f''''(ubound) = f^(2*beta)(ubound) = ucond
%   The convergence as beta increases of the maternqr methods can then
%   be compared between f, which will not in general satisfy the boundary
%   conditions, and g, which will. As the only difference between f and g
%   is the terms added by this function, this may give some insight into
%   how satisfaction of the boundary conditions affects maternqr
%   convergence.
syms x
if beta == -1
    g = f; % no boundary conditions will be forced at all; just return f
    return
elseif beta == 0
    d = f; % just compensating for the deviation from boundary values of f
           % itself
else
    d = diff(f,2*beta); % symbolically find the highest order derivative
end

%--------------------------------------------------------------------------
% Construct "adjustment" line that will correct any deviation from desired 
% boundary values in the highest derivative
yone = subs(d-lbval,x,sym(lbound)); % find coordinates of the two
xone = sym(lbound);                  % points through which our adjustment line
ytwo = subs(d-ubval,x,sym(ubound));  % will pass
xtwo = sym(ubound);

m = (ytwo-yone)/(xtwo-xone); % find the slope of the adjustment line

% Construct (polynomial-vector representation of) adjustment line:
adjustmentTerm = sym2poly( m *(x - xone) + yone );
%--------------------------------------------------------------------------

% symbollically integrate the adjustment line 2*beta times to find the term
% we must add to f so that its highest derivative satisfyies the boundary
% conditions
for i = 1:2*beta
    adjustmentTerm = polyint(adjustmentTerm);
end

newf = f - poly2sym(adjustmentTerm);  % this new function's (2*beta)th derivative
                            % satisfies the boundary conditions
                            
if beta == 0
    g = symfun(newf,x);
else % recursive call on the new function with beta decremented
    g = genBoundarySatisfyingFunction(newf, beta-1, lbound, ubound, lbval, ubval);
end

end


function [ g ] = forceBCsatADD( f, lbsat, rbsat, lbound, rbound, lbval, rbval )
%GENBOUNDARYSATISFYINGFUNCTION
% Takes a symbolic f(x) and generates the function g(x) = f(x)+(terms
%that cause f and all even derivatives of f up to the (2*beta)th derivative
%to equal 0 at the at the boundary of [lbound,ubound])
%   Returns a function g based on f that satisfies the boundary 
%   conditions
%      f(lbound) = f''(lbound) = f''''(lbound) = f^(2*beta)(lbound) = lcond
%      f(ubound) = f''(ubound) = f''''(ubound) = f^(2*beta)(ubound) = ucond
%   The convergence as beta (not N!) increases of the maternqr methods can 
%   be compared between f, which will not in general satisfy the boundary
%   conditions, and g, which will. As the only difference between f and g
%   is the terms added by this function, this may give some insight into
%   how satisfaction of the boundary conditions affects maternqr
%   convergence.
syms x
if lbsat == rbsat
    curSatDegree = lbsat;
    curlbval = lbval; % Force both sides to satisfy boundary conds
    currbval = rbval;
    newlbsat = lbsat - 1; % decrement both satisfaction degrees
    newrbsat = rbsat - 1;
    if lbsat == -1 && rbsat == -1
        g = f;  % no boundary conditions will be forced at all; 
                % just return f
    return
    elseif lbsat == 0 && rbsat == 0
        d = f;  % just compensating for the deviation from boundary values 
                % of f itself
        
    else % symbolically find the highest order derivative:
        d = diff(f,2*lbsat); 
    end
elseif lbsat > rbsat
    d = diff(f,2*lbsat); % symbolically find the highest order derivative
    curlbval = lbval; % DO force highest-derivative bval satistfaction on 
                      % LEFT boundary
    currbval = d; % DON'T force highest-derivative bval satistfaction on 
                  % RIGHT boundary
    newlbsat = lbsat - 1; % decrement only left satisfaction degree, 
    newrbsat = rbsat;     % since that's all we've forced
    curSatDegree = lbsat;
elseif lbsat < rbsat
    d = diff(f,2*rbsat); % symbolically find the highest order derivative
    curlbval = d; % DON'T force highest-derivative bval satistfaction on 
                  % LEFT boundary
    currbval = rbval; % DO force highest-derivative bval satistfaction on 
                      % RIGHT boundary
	newlbsat = lbsat;     % decrement only right satisfaction degree, 
    newrbsat = rbsat - 1; % since that's all we've forced
    curSatDegree = rbsat;
end

%--------------------------------------------------------------------------
% Construct "adjustment" line that will correct any deviation from desired 
% boundary values in the highest derivative
yleft = subs(d-curlbval,x,sym(lbound));  % find coordinates of the two
xleft = sym(lbound);                     % points through which our adjustment
yright = subs(d-currbval,x,sym(rbound)); % line will pass
xright = sym(rbound);

m = (yright-yleft)/(xright-xleft); % find the slope of the adjustment line

% Construct (polynomial-vector representation of) adjustment line:
adjustmentTerm = sym2poly( m *(x - xleft) + yleft );
%--------------------------------------------------------------------------

% symbollically integrate the adjustment line 2*beta times to find the term
% we must add to f so that its highest derivative satisfyies the boundary
% conditions
for i = 1:2*curSatDegree
    adjustmentTerm = polyint(adjustmentTerm);
end

% this new function's (2*beta)th derivative
% satisfies the boundary conditions
newf = f - poly2sym(adjustmentTerm); 
                            
if curSatDegree == 0
    g = symfun(newf,x);
else % recursive call on the new function with beta decremented
    g = forceBCsatADD(newf, newlbsat, newrbsat, lbound, rbound, lbval, rbval);
end

end


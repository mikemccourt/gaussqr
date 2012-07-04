function g = forceBCsatMULT( f, lbsat, rbsat, lbound, rbound)
%FORCEBCSATMULT
% Takes a symbolic formula for a function f, and returns a symbolic
% function g that satisfies the boundary conditions
%       (ith derivative of g)(lbound) = 0
%       (jth derivative of g)(rbound) = 0
% for 0 <= i <= 2*lbsat and 0 <= j <= 2*rbsat. This is done by multiplying
% by polynomials with zeros at lbound and rbound.
    syms x
    g = symfun(f * (x - lbound)^(2*lbsat+1) * (x - rbound)^(2*rbsat+1),x);
end
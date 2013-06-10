function f = sinc(x)
% function f = sinc(x)
% This function evaluates the function sinc(x) = sin(x)/x
% It defines sinc(0)=1
%
% Inputs: x - data locations
% Outputs: f - sinc(x) values
%
% NOTE: As with any sin(x) evaluation, extreme values of x may produce
% inaccurate results, though this is dampened by the 1/x term.

f = ones(size(x));
nz = find(x~=0);
f(nz) = sin(pi*x(nz))./(pi*x(nz));

% f = sinc(x)
%
% Defines sinc function

function f = sinc(x)

f = ones(size(x));
nz = find(x~=0);
f(nz) = sin(pi*x(nz))./(pi*x(nz));

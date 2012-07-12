function out = testkernel(x)
% TESTKERNEL
% A 1-D slice of the beta=2,epsilon>0 compact Matern kernel.
% This kernel has a discontinuity in the 3rd derivative at x = y.

% set the plane with which to slice our 2-D kernel to obtain a 1-D curve:
y = .5; % this is where the discontinuity in the 3rd derivative will occur
eps = 1;
L = 1;

out = zeros(length(x),1);
for i = 1:length(x)
    if (0 <= x(i) )&& (x(i) < y)
        out(i) = (1/(2*eps^3*(-1+cosh(2*L*eps)+sinh(2*L*eps))^2))*(4*x(i)*eps*cosh(x(i)*eps)*sinh(L*eps)*(cosh(2*L*eps)+sinh(2*L*eps))*sinh((L-y)*eps)+sinh(x(i)*eps)*(cosh(y*eps)-sinh(y*eps))*((1-2*L*eps+y*eps)*(cosh(2*L*eps)+sinh(2*L*eps))-(1+y*eps)*(cosh(4*L*eps)+sinh(4*L*eps))+(-1+y*eps)*(cosh(2*y*eps)+sinh(2*y*eps))+(1+2*L*eps-y*eps)*(cosh(2*(L+y)*eps)+sinh(2*(L+y)*eps))));
    elseif ((y <= x(i)) && (x(i) <= 1))
        out(i) = (1/(2*eps^3*(-1+cosh(2*L*eps)+sinh(2*L*eps))^2))*(4*y*eps*cosh(y*eps)*sinh(L*eps)*(cosh(2*L*eps)+sinh(2*L*eps))*sinh((L-x(i))*eps)+sinh(y*eps)*(cosh(x(i)*eps)-sinh(x(i)*eps))*((1-2*L*eps+x(i)*eps)*(cosh(2*L*eps)+sinh(2*L*eps))-(1+x(i)*eps)*(cosh(4*L*eps)+sinh(4*L*eps))+(-1+x(i)*eps)*(cosh(2*x(i)*eps)+sinh(2*x(i)*eps))+(1+2*L*eps-x(i)*eps)*(cosh(2*(L+x(i))*eps)+sinh(2*(L+x(i))*eps))));
    else
        error('An x value is not in the range (0,1)');
    end
end

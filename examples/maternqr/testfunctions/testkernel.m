function out = testkernel(x)
% TESTKERNEL
% A 1-D slice of the beta=2,epsilon>0 compact Matern kernel.
% This kernel has a discontinuity in the 3rd derivative at x = y.

% set the plane with which to slice our 2-D kernel to obtain a 1-D curve:
y = .5; % this is where the discontinuity in the 3rd derivative will occur

out = zeros(length(x),1);
for i = 1:length(x)
    if (0 <= x(i) )&& (x(i) < y)
        out(i) = (1/(2*(exp(2)-1))^2)*((exp(2)-1)*x(i)*cosh(x(i))*((exp(2)-1)*cosh(y)-(1+exp(2))*sinh(y))+sinh(x(i))*(cosh(y)-sinh(y))*(exp(2)*(y-1)-exp(4)*(y+1)+(y-1)*(sinh(2*y)+cosh(2*y))-(y-3)*sinh(2*(y+1))+cosh(2*(y+1))));
    elseif ((y <= x(i)) && (x(i) <= 1))
        out(i) = (1/(2*(exp(2)-1))^2)*((exp(2)-1)*y*cosh(y)*((exp(2)-1)*cosh(x(i))-(1+exp(2))*sinh(x(i)))+sinh(y)*(cosh(x(i))-sinh(x(i)))*(exp(2)*(x(i)-1)-exp(4)*(x(i)+1)+(x(i)-1)*(sinh(2*x(i))+cosh(2*x(i)))-(x(i)-3)*sinh(2*(x(i)+1))+cosh(2*(x(i)+1))));
    else
        error('An x value is not in the range [0,1]');
    end
end

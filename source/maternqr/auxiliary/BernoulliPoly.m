function b = BernoulliPoly(n,x)
% function b = BernoulliPoly(n,x)
% This function evaluates the Bernoulli polynomial of order n at the vector
% specified in x.  This is NOT a smart or well-tested version, I just wrote
% it up to run some quick tests for small n.
%
% Inputs: n - order of the polynomial
%         x - column vector of locations at which to evaluate
%
% Outputs: b - the Bernoulli polynomial of order n at x

N = size(x,1);

b = zeros(N,1);
for m=0:n
    sp = repmat((-1).^(0:m).*exp(gammaln(m+1)-gammaln((0:m)+1)-gammaln(m-(0:m)+1)),N,1).*...
        (repmat(x,1,m+1)+repmat(0:m,N,1)).^n;
    b = b + sum(sp,2)/(m+1);
end
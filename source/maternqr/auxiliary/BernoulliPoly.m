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

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
bn = GAUSSQR_PARAMETERS.BERNOULLI_NUMBERS;

if(n>length(bn)-1)
    error('Bernoulli polynomials of degree>%d must be computed with toolbox',length(bn)-1)
end

N = size(x,1);

B = repmat(bn(1:n+1).*exp(gammaln(n+1)-gammaln((0:n)+1)-gammaln(n-(0:n)+1)),N,1);
X = repmat(x,1,n+1).^repmat(n-(0:n),N,1);
b = sum(B.*X,2);
function gam = eulerGamma(NDig)
% eulerGamma: Compute the Euler-Mascheroni constant (gamma)
% usage: gam = eulerGamma(NDig)
%
% Uses the Bessel function method from:
% http://numbers.computation.free.fr/Constants/Gamma/gamma.pdf
%
% See also:
% http://en.wikipedia.org/wiki/Euler?Mascheroni_constant
%
% arguments: (input)
%  NDig - the number of correct digits in the value of gamma
%
% arguments: (output)
%  gam - the Euler-Mascheroni constant (gamma) as an HPF number
%
%

% check the argument(s)
if nargin == 0
  NDig = DefaultNumberOfDigits;
elseif nargin > 1
  error('EULERGAMMA:improperarguments', ...
    'Zero or one arguments only')
elseif ~ismember(numel(NDig),[1,2]) || any(NDig ~= round(NDig))
  error('EULERGAMMA:improperarguments', ...
    'NDig must be scalar integer or vector of length 2')
end

% get alpha such that: alpha*(log(alpha) - 1) == 1
alpha = fzero(@(alpha) alpha*(log(alpha) - 1) - 1,[3,5]);

% Choose a value of n, such that the error is on the order of eps
% for the given total number of digits. For the scheme used here,
% the error is O(exp(-4*n)). Therefore,
n = ceil(sum(NDig)./(4*log10(exp(1))));

% The nth harmonic number. H0 = 0
Hn = hpf(0,NDig);
An = hpf(0,NDig);
Bn = hpf(1,NDig);
C = hpf(1,NDig);
n_hpf = hpf(n,NDig);
for k = 1:floor(alpha*n)
  Hn = Hn + reciprocal(hpf(k,NDig));
  C = C.*(n_hpf./k).^2;
  Bn = Bn + C;
  An = An + C*Hn;
end

gam = An./Bn - log(hpf(n,NDig));

end




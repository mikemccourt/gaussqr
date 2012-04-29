function [Happx,negInd] = HermiteAppx(n,x,logopt)
% This function approximates Hermite polynomials for large n.
% INPUTS:
%   n - Order of the Hermite polynomial
%   x - Column vector of data points
%   logopt - Set to 1 to return the log of the polynomial
% OUTPUTS:
%   Happx - The Hermite Approximation
%   negInd - Indices which are actually negative
%            Only used if logopt=1 to avoid complex logarithms

[N,c] = size(x);
if c~=1
    error(sprintf('x must be a column vector: mr=%d, mc=%d',N,c))
end
if not(exist('logopt'))
    logopt=0;
elseif not(logopt==0 || logopt==1)
    logopt=0;
end
Happx = zeros(size(x));

absx = abs(x);
s2n = sqrt(2*n);
transWidth = .03; % Maybe this should be a variable width
innerInd = find(absx<s2n*(1-transWidth))';
outerInd = find(absx>s2n*(1+transWidth))';
transInd = setdiff(1:N,[innerInd,outerInd]);

[innerAppx,innernegLog] = HermiteAppxInner(n,absx(innerInd),logopt);
[transAppx,transnegLog] = HermiteAppxTrans(n,absx(transInd),logopt);
outerAppx = HermiteAppxOuter(n,absx(outerInd),logopt);

Happx([innerInd,transInd,outerInd]) = [innerAppx;transAppx;outerAppx];
negInd = zeros(size(x));
if (logopt==0)
    if (mod(n,2)==1)
        negInd = find(x<0);
        Happx(negInd) = -Happx(negInd);
    end
else
    if mod(n,2)==1 % Odd power means odd polynomial
        negInd = x<0;
    end
    negInd([innerInd,transInd]) = xor(negInd([innerInd,transInd]),[innernegLog;transnegLog]);
end

end

% Note below that sqrt(2/cos(t)) = sqrt(2/sqrt(1-x^2/2n)) -> 1/2log(2)-1/4log(1-x^2/2n)
function [Happx,negLog] = HermiteAppxInner(n,x,logopt) % Only valid in abs(x)<sqrt(2n)
xos2n = x/sqrt(2*n);
xos2nS = xos2n.^2;
omxos2nS = 1 - xos2nS;
t = asin(xos2n);
if logopt==1
    costerm = n*(xos2n.*sqrt(omxos2nS)+t-pi/2)+t/2;
    negLog = mod(costerm-pi/2,2*pi)<pi;
    Happx = 1/2*log(2)-1/4*log(omxos2nS) + ...
            n/2*(log(2*n)+1-2*omxos2nS) + ...
            log(abs(cos(costerm))+eps);
else
    negLog = zeros(size(x));
    Happx = sqrt(2)*(omxos2nS).^(-1/4).* ... 
            exp(n/2*(log(2*n)+1-2*omxos2nS)).* ...
            cos(n*(xos2n.*sqrt(omxos2nS)+t-pi/2)+t/2);
end
end

% This has limited applicability for extremely large values of x
% The asymptotic form is airy(x) -> exp(-2/3x^(3/2))/(2*sqrt(pi)*x^.25)
% I haven't implemented it yet, but I may
%
% Note: For real x<0, airy(x) may have a trivial complex component
%       That value is disregarded as it should not exist
%
% Weird: The airy function acts differently on an empty matrix
%        than exp does.  This is why there is that weird start
function [Happx,negLog] = HermiteAppxTrans(n,x,logopt) % Only valid near x=sqrt(2n)
if length(x)==0
    Happx = [];
    negLog = [];
else
    if logopt==1
        airyterm = real(airy(sqrt(2)*n^(1/6)*(x-sqrt(2*n))));
        negLog = airyterm<0;
        Happx = 1/2*log(2*pi)+1/6*log(n)+n/2*log(2*n)-3/2*n + ...
            sqrt(2*n)*x + ...
            log(abs(airyterm)+eps);
    else
        negLog = zeros(size(x));
        Happx = sqrt(2*pi)*n^(1/6)*exp(n/2*log(2*n)-3/2*n)* ...
            exp(sqrt(2*n)*x).* ...
            real(airy(sqrt(2)*n^(1/6)*(x-sqrt(2*n))));
    end
end
end

function Happx = HermiteAppxOuter(n,x,logopt) % Only valid in x>sqrt(2n)
s = sqrt(x.^2-2*n);
if logopt==1
    Happx = 1/2*(-log(2)+log(1+x./s)+(x.^2-s.*x-n)) + n*log(s+x);
else
    Happx = sqrt(1/2*(1+x./s)).*exp(1/2*(x.^2-s.*x-n)+n*log(s+x));
end
end
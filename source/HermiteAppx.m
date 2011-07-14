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
    error('Only a column vector can be passed to HermiteAppx')
end
if not(exist('logopt'))
    logopt=0;
elseif not(logopt==0 || logopt==1)
    logopt=0;
end

t = asin(x/sqrt(2*n));
Happx = sqrt(2./cos(t)).*exp(n/2*(log(2*n)-cos(2*t))).*cos(n*(sin(2*t)/2+t-pi/2)+t/2);
cosTerm = cos(n*(sin(2*t)/2+t-pi/2)+t/2);
negLog = find(cosTerm<0);
logHappx = 1/2*log(2)-1/4*(log(1-x.^2/(2*n)))+log(abs(cosTerm))+n/2*(log(2*n)-cos(2*t));
Happx = exp(logHappx);
Happx(negLog)=-Happx(negLog);

if mod(n,2)==1 % Odd power means odd polynomial
    negInd = find(x<0);
else
    negInd = [];
end
absx = abs(x);
s2n = sqrt(2*n);
transWidth = .1; % Maybe this should be a variable width
innerInd = find(absx<s2n*(1-transWidth));
outerInd = find(absx>s2n*(1+transWidth));
transInd = setdiff(1:N,[innerInd,outerInd]);

if length(innerInd>0)
    [innerAppx,innernegLog] = HermiteAppxInner(n,x(innerInd),logopt);
end
if length(outerInd>0)
    [outerAppx,outernegLog] = HermiteAppxOuter(n,x(outerInd),logopt);
end
if length(transInd>0)
    [transAppx,transnegLog] = HermiteAppxTrans(n,x(transInd),logopt);
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
    modcosterm = mod(costerm-pi,2*pi);
    negLog = (modcosterm<=-pi/2) + (modcosterm>=pi/2);
    negLog = 1 - negLog;
    Happx = 1/2*log(2)-1/4*log(omxos2nS) + ...
            n/2*(log(2*n)+1-2*omxos2nS) + ...
            log(abs(cos(costerm)));
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
function [Happx,negLog] = HermiteAppxTrans(n,x,logopt) % Only valid near abs(x)=sqrt(2n)
if logopt==1
    airyterm = real(airy(sqrt(2)*n^(1/6)*(x-sqrt(2*n))));
    negLog = airyterm<0;
    Happx = 1/2*log(2*pi)+1/6*log(n)+n/2*log(2*n)-3/2*n + ...
            sqrt(2*n)*x + ...
            log(abs(airyterm));
else
    negLog = zeros(size(x));
    Happx = sqrt(2*pi)*n^(1/6)*exp(n/2*log(2*n)-3/2*n) ...
            exp(sqrt(2*n)*x).* ...
            airy(sqrt(2)*n^(1/6)*(x-sqrt(2*n)));
end
end

function [Happx,negLog] = HermiteAppxOuter(n,x,logopt) % Only valid in abs(x)>sqrt(2n)
negLog = [];
t = asin(x/sqrt(2*n));
Happx = sqrt(2./cos(t)).*exp(n/2*(log(2*n)-cos(2*t))).*cos(n*(sin(2*t)/2+t-pi/2)+t/2);
end
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
transWidth = .1;
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

function [Happx,negLog] = HermiteAppxInner(n,x,logopt) % Only valid in abs(x)<sqrt(2n)
negLog = [];
t = asin(x/sqrt(2*n));
Happx = sqrt(2./cos(t)).*exp(n/2*(log(2*n)-cos(2*t))).*cos(n*(sin(2*t)/2+t-pi/2)+t/2);
end

function [Happx,negLog] = HermiteAppxOuter(n,x,logopt) % Only valid in abs(x)>sqrt(2n)
negLog = [];
t = asin(x/sqrt(2*n));
Happx = sqrt(2./cos(t)).*exp(n/2*(log(2*n)-cos(2*t))).*cos(n*(sin(2*t)/2+t-pi/2)+t/2);
end

function [Happx,negLog] = HermiteAppxTrans(n,x,logopt) % Only valid near abs(x)=sqrt(2n)
negLog = [];
t = asin(x/sqrt(2*n));
Happx = sqrt(2./cos(t)).*exp(n/2*(log(2*n)-cos(2*t))).*cos(n*(sin(2*t)/2+t-pi/2)+t/2);
end
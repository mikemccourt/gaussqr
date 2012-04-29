function Marr = gqr_formMarr(M,Mmax,maxlength)
% function Marr = gqr_formMarr(M,Mmax,maxlength)
%
% Inputs:
%   M - an integer vector that tells the max power in each dimension.
%       If you pass M=[0;...;0] this will set each M to Mmax
%   Mmax - an optional integer giving the max power that epsilon should see.
%          This may be omitted by passing []
%   maxlength - optional, and limits size(Marr,2)<maxlength
%               This prevents Marr from being too long

M = M(:)'; % Make it a row vector to make my life easier
d = length(M);
if sum(M)==0
    if not(exist('Mmax'))
        Mmax = 1;
    elseif length(Mmax)==0
        if not(exist('maxlength'))
            Mmax = 1;
        else
            if maxlength<inf
                Mmax = maxlength;
            else
                error('Incompatible Mmax=%g and maxlength=%g values',Mmax,maxlength)
            end
        end
    end
    M = ones(1,d)*Mmax;
end
if not(exist('Mmax'))
    Mmax = prod(M);
elseif length(Mmax)==0
    Mmax = prod(M);
end
if not(exist('maxlength'))
    maxlength = inf;
end

Marr = zeros(d,1);
x = [1];
Msum = 1;
while Msum<=Mmax && size(Marr,2)<maxlength
    for k=x
        z = repmat(Marr(:,k),1,d)+eye(d);
        Marr = [Marr,z(:,find(sum(z)<M))];
    end
    [b,m,n] = unique(Marr','rows');
    Marr = Marr(:,sort(m));
    x = find(sum(Marr,1)==Msum);
    Msum = Msum + 1 + Mmax*(length(x)==0);
end
if size(Marr,2)>maxlength
    Marr = Marr(:,1:maxlength);
end
Marr = Marr + 1; % For the off-by-one to be corrected
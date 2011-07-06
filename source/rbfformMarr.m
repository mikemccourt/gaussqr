function Marr = rbfformMarr(M,Mmax)
% function Marr = rbfformMarr(M,Mmax)
% M should be a vector that tells the max length in each dimension.
% Mmax is an optional input giving the max power that epsilon should see.
% If you pass this M=[0;...;0] this will set each M to Mmax
M = M(:)'; % Make it a row vector to make my life easier
d = length(M);
if sum(M)==0
    if not(exist('Mmax'))
        Mmax = 1;
    end
    M = ones(1,d)*Mmax;
end
if not(exist('Mmax'))
    Mmax = prod(M);
end
Marr = zeros(d,1);
x = [1];
Msum = 1;
while Msum<=Mmax
    for k=x
        z = repmat(Marr(:,k),1,d)+eye(d);
        Marr = [Marr,z(:,find(sum(z)<M))];
    end
    [b,m,n] = unique(Marr','rows');
    Marr = Marr(:,sort(m));
    x = find(sum(Marr,1)==Msum);
    Msum = Msum + 1 + Mmax*(length(x)==0);
end
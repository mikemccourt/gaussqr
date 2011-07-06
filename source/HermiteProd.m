function Hx = HermiteProd(m,x,logoption)
% This function takes a vector m and a matrix x
% [mr mc] = size(m)
% [xr xc] = size(x)
% We need xc==mr since that is the number of dimensions
% We need mc==1 since we can only evaluate one product at a time
% xr is the number of points we need to evaluate at
% If logoption is activated, the value returned will be log of Hx
%   Note that with logoption on you will have complex Hx
%   You should see real(exp(HermiteProd(m,x,1)))=HermiteProd(m,x) to near machine precision 

M = {[1] [2 0] [4 0 -2] [8 0 -12 0] [16 0 -48 0 12] [32 0 -160 0 120 0] [64 0 -480 0 720 0 -120] ...
     [128 0 -1344 0 3360 0 -1680 0] [256 0 -3584 0 13440 0 -13440 0 1680] [512 0 -9216 0 48384 0 -80640 0 30240 0] ...
     [1024 0 -23040 0 161280 0 -403200 0 302400 0 -30240]};

[mr mc] = size(m);
[xr xc] = size(x);

if xc~=mr
    error(sprintf('Dimension mismatch: xc=%d, mr=%d',xc,mr))
end
if mc~=1
    error(sprintf('m must be a column vector: mr=%d, mc=%d',mr,mc))
end
if nargin==3 && (logoption~=0 || logoption~=false)
    logoption=1;
else
    logoption=0;
end

for k=1:length(m)
    if m(k)<=length(M)-1
        H = polyval(M{m(k)+1},x(:,k));
    else
        H = polyval(HermitePoly(m(k)),x(:,k));
    end
    if logoption==1
        if k==1
            Hx = log(H+eps); % For log(0)
        else
            Hx = Hx + log(H+eps); % For log(0)
        end
    else
        if k==1
            Hx = H;
        else
            Hx = Hx.*H;
        end
    end
end
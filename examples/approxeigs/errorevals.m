function [errr,errc] = errorevals(D,C)
    realr = @(n) 1./((n.^2)*pi^2);
    s = size(D);
    n = [1: s(1)]';
    Rer = realr(n);
    errr = abs( D - Rer)./Rer;
    errc = C;%hack here
    
    %little confused about this part. It seems errcompute cannot compute
    %the error between two fucntion. Should we choose some pointes to
    %figure the error?
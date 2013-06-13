%function [errr,errc] = errorevals(D,C,M)
%Evaluates the error for the eigenvalue and eigenfunction
%
%Inputs : D - a vector of sorted eigenvalue
%         C - a matrix of eigenfunction 
%         M - different ways to measure the error
%Output : errr - a vector of eigenvalues' relative error
%       : errc - a vector of eigenfunctions' relative error
function [errr,errc] = errorevals(D,C,M)
    % compute eigenvalues' relative error 
    realr = @(n) 1./((n.^2)*pi^2);
    s = size(D);
    n = [1: s(1)]';
    Rer = realr(n);
    errr = abs( D - Rer)./Rer;
    % compute eigenfunctions' realtive error
    sc = size(C);
    j = sc(2); %how many columns
    rl = sc(1); %how many rows 
    errc = zero(rl,1);
    switch M
        case 1 %abs(normone(realfunction)-normone(apporximatefunction))
            for i =1: j
                CV = C(:,i);
                CV = CV*(CV(2)/abs(CV(2)));
                k = (0 : rl-1);
                f = @(x)(x.^k)*CV;
                maxf = max(f(0 : 0.002 : 1));
                ff = @(x) abs((x.^k) * CV / maxf * sqrt(2));
                errc(i) = abs( 2*sqrt(2)/pi - quad(ff,0,1));
            end
        case 2 %normone(realfunction - apporximatefunction)
            for i =1: j
                CV = C(:,i);
                CV = CV*(CV(2)/abs(CV(2)));
                k = (0 : rl-1);
                f = @(x) (x.^k)*CV;
                maxf = max(f(0 : 0.002 : 1));
                ff = @(x) abs(sqrt(2)*sin(pi*i*x) - (x.^k) * CV / maxf * sqrt(2));
                errc(i) = quad(ff,0,1);
            end
        case 3 %boundary value
             for i =1: j
                CV = C(:,i);
                CV = CV*(CV(2)/abs(CV(2)));
                k = (0 : rl-1);
                f = @(x)(x.^k)*CV;
                maxf = max(f(0 : 0.002 : 1));
                ff = @(x) (x.^k) * CV / maxf * sqrt(2);
                errc(i) = abs(ff(0))/sqrt(2); %max value of function ff should be sqrt(2)
             end
        case 4 %collection points
             for i =1: j
                CV = C(:,i);
                CV = CV*(CV(2)/abs(CV(2)));
                k = (0 : rl-1);
                f = @(x)(x.^k)*CV;
                maxf = max(f(0 : 0.002 : 1));
                ff = @(x) abs(sqrt(2)*sin(i*pi*x)-(x.^k) * CV / maxf * sqrt(2));
                errc(i) =sum(ff(0:.002:1))/500 ; %sommthing wrong here
             end
            
        otherwise errc = errc;
    end
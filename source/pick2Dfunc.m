function [yf,fstr] = pick2Dfunc(fopt)
% function [x,spacestr] = pick2Dpoints(a,b,N,spaceopt,ep)
% Must pass a=[D1 lower bound,   D2 lower bound]
%           b=[D1 upper bound,   D2 upper bound]
%           N=[num points in D1, num points in D2]
% Note that ep is only needed for 'inner' points
% This assumes that a, b, are 2 element row vectors
% N should be a 2-vector, but if N is scalar, we reset N=[N,N]
% Note that the Halton points only really need prod(N)

switch lower(fopt)
    case 'poly'
        fstr = 'f(x,y)=x+y/2+xy/4';
%         yf = @(x) (1-x(:,1).^2).*(1-x(:,2).^2);
        yf = @(x) x(:,1)/1+x(:,2)/2 + x(:,1).*x(:,2)/4;
    case 'sin'
        fstr = 'f(x,y)=cos((x^2+y^2))';
%         yf = @(x) cos((x(:,1).^2+x(:,2).^2)/100);
        yf = @(x) cos((x(:,1)+x(:,2)));
    case 'runge'
        fstr = 'f(x,y)=1/(1+x^2+y^2)';
        yf = @(x) 1./(1+(x(:,1).^2+x(:,2).^2));
    case 'tanh'
        fstr = 'f(x,y)=tanh(9(y-x)+1)/(tanh(9)+1)';
        yf = @(x) (tanh(9*(x(:,2)-x(:,1)))+1)/(tanh(9)+1);
        fstr = 'f(x,y)=tanh(x+y)';
        yf = @(x) tanh(x(:,1)+x(:,2));
    case 'franke'
        fstr = 'f(x,y) = Frankes function';
        yf = @(x) franke(x(:,1)/1,x(:,2)/1);
end
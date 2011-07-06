function [yf,fstr] = pickfunc(fopt,dim)
% function [yf,fstr] = pickfunc(fopt,dim)

switch dim
    case 1
        switch lower(fopt)
            case 'poly'
                fstr = 'y(x)=x^3-2x^2+.5x';
                yf = @(x) x.^3-2*x.^2+.5*x;
            case 'sin'
                fstr = 'y(x)=sin(x/2)-2cos(x)+4sin(\pi x)';
                yf = @(x) sin(x/2)-2*cos(x)+4*sin(pi*x);
            case 'runge'
                fstr = 'y(x)=1/(1+x^2)';
                yf = @(x) 1./(1+x.^2);
            case 'sinh'
                fstr = 'y(x)=sinh(x)/(1+cosh(x))';
                yf = @(x) sinh(x)./(1+cosh(x));
        end
    case 2
        switch lower(fopt)
            case 'poly'
                fstr = 'f(x,y)=x+y/2+xy/4';
                yf = @(x) x(:,1)/1+x(:,2)/2 + x(:,1).*x(:,2)/4;
            case 'sin'
                fstr = 'f(x,y)=cos((x^2+y^2))';
                yf = @(x) cos((x(:,1)+x(:,2)));
            case 'runge'
                fstr = 'f(x,y)=1/(1+x^2+y^2)';
                yf = @(x) 1./(1+(x(:,1).^2+x(:,2).^2));
            case 'tanh'
                fstr = 'f(x,y)=tanh(x+y)';
                yf = @(x) tanh(x(:,1)+x(:,2));
            case 'franke'
                fstr = 'f(x,y) = Frankes function';
                yf = @(x) franke(x(:,1)/1,x(:,2)/1);
        end
    otherwise
        error('Can only consider 1D and 2D functions')
end
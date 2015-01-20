function [yf,fstr] = pickfunc(fopt,dim)
% function [yf,fstr] = pickfunc(fopt,dim)
% For the optics functions in 2D, only consider the domain [-1,1] to get
% the real problem because the scaling is accounted for within this
% function.
%
% PROGRAMMERS NOTE: It may be better to eliminate/alter the dimension
% option, as I have tried to do for some of the surrogate models.  Part of
% what this would require is a dictionary of acceptable functions.

% Prevent any uppercase/lowercase issues from popping up
fopt = lower(fopt);

% Identify if there are special functions
if any(strcmp(fopt,{'piston','piston_scaled'}))
    dim = 7;
end

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
            case 'exp'
                fstr = 'y(x)=10e^{-x^2}+x^2';
                yf = @(x) 10*exp(-x.^2)+x.^2;
            case 'cosexp'
                fstr = 'y(x)=cos(x)+e^{-(x-1)^2}-e^{-(x+1)^2}';
                yf = @(x) cos(x)+exp(-(x-1).^2)-exp(-(x+1).^2);
            case 'log'
                fstr = 'y(x)=xlog(1+x^2)';
                yf = @(x) x.*log(1+x.^2);
            case 'qian'
                fstr = 'y(x)=exp((x-1)^2) sin(exp((x-1)^2))';
                yf = @(x) exp((x-1).^2).*sin(exp((x-1).^2));
            otherwise
                error('No such function %s exists in %d dimensions',fopt,dim)
        end
    case 2
        switch lower(fopt)
            case 'poly'
                fstr = 'f(x,y)=x+y/2+xy/4';
                yf = @(x) x(:,1)/1+x(:,2)/2 + x(:,1).*x(:,2)/4;
            case 'sin'
                fstr = 'f(x,y)=cos((x+y))';
                yf = @(x) cos((x(:,1)+x(:,2)));
            case 'runge'
                fstr = 'f(x,y)=1/(1+x^2+y^2)';
                yf = @(x) 1./(1+(x(:,1).^2+x(:,2).^2));
            case 'tanh'
                fstr = 'f(x,y)=tanh(x+y)';
                yf = @(x) tanh(x(:,1).^2+x(:,2).^2);
            case 'franke'
                fstr = 'f(x,y) = Frankes function';
                yf = @(x) franke(x(:,1),x(:,2));
            case 'kxy' % for [-5,5]^2
                fstr = 'f(x,y) = KXY';
                yf = @(x) kxy(5*x);
            case 'ksa1' % for [-1,1]^2
                fstr = 'f(x,y) = KSA1';
                yf = @(x) ksa1(x);
            case 'ksa2' % for disk with radius 20
                fstr = 'f(x,y) = KSA2';
                yf = @(x) ksa2(20*x);
            otherwise
                error('No such function %s exists in %d dimensions',fopt,dim)
        end
    case 7
        switch(fopt)
            case 'piston'
                fstr = 'f(x7) = piston motion';
                yf = @(x) piston(x(:,1),x(:,2),x(:,3),x(:,4),x(:,5),x(:,6),x(:,7));
            case 'piston_scaled'
                fstr = 'f(x7) = piston motion [0,1]^7';
                lb = [30,.005,.002,1000,90000 ,290,340];
                ub = [60,.020,.010,5000,110000,296,360];
                yf = @(x) piston(x(:,1)*(ub(1)-lb(1)) + lb(1),...
                                 x(:,2)*(ub(2)-lb(2)) + lb(2),...
                                 x(:,3)*(ub(3)-lb(3)) + lb(3),...
                                 x(:,4)*(ub(4)-lb(4)) + lb(4),...
                                 x(:,5)*(ub(5)-lb(5)) + lb(5),...
                                 x(:,6)*(ub(6)-lb(6)) + lb(6),...
                                 x(:,7)*(ub(7)-lb(7)) + lb(7));
            otherwise
                error('No such function %s exists in %d dimensions',fopt,dim)
        end
    otherwise
        error('No functions appear in %d dimensions',dim)
end
end

function y = kxy(x)
% used on [-5,5]^2 in [Jester/Menke/Urban (2011)]
    rx = -2.15e-1;
    ry = 1.22e-1;
    kx = -1;
    ky = -1; 
    xx = x(:,1);
    yy = x(:,2);
    y = (rx*xx.^2+ry*yy.^2)./(1+sqrt(1-(1+kx)*rx^2*xx.^2-(1+ky)*ry^2*yy.^2)) ...
        - 4.05e-4*xx.^4 - 8.13e-4*xx.^2.*yy.^2 + 5.73e-4*yy.^4 ...
        - 4.59e-6*xx.^6 + 1.14e-5*xx.^4.*yy.^2 + 9.64e-6*xx.^2.*yy.^4 + 4.45e-7*yy.^6 ...
        - 2.69e-9*xx.^8 - 7.96e-8*xx.^6.*yy.^2 - 8.79e-8*xx.^4.*yy.^4 - 9.16e-8*xx.^2.*yy.^6 + 2.43e-8*yy.^8;
end
function y = ksa1(x)
% used on [-1,1]^2 in [Jester/Menke/Urban (2011)]
    r2 = x(:,1).^2+x(:,2).^2;
    rho = 1;
    kappa = -1;
    y = rho*r2./(1+sqrt(1-(1+kappa)*rho^2*r2)) + 1.5*r2.^2 ...
            - 7e-1*r2.^3 + 5e-1*r2.^4 - 5e-2*r2.^5;
end
function y = ksa2(x)
% used on disk of radius 20 in [Jester/Menke/Urban (2011)]
% here computes on a disk of radius 1
    r2 = x(:,1).^2+x(:,2).^2;
    rho = -3.87e-2;
    kappa = 0;
    y = rho*r2./(1+sqrt(1-(1+kappa)*rho^2*r2)) - 4.17e-6*r2.^2 + 4.71e-9*r2.^3 ...
        - 4.94e-12*r2.^4 - 5.42e-15*r2.^5 - 4.98e-18*r2.^6 - 1.22e-20*r2.^7;
end

function f = franke(x,y)

f = .75*exp(-.25*((9*x-2).^2+(9*y-2).^2))+ ...
                     .75*exp(-1/49*(9*x+1).^2-.1*(9*y+1).^2)+ ...
                     .5*exp(-.25*((9*x-7).^2+(9*y-3).^2))- ...
                     .2*exp(-(9*x-4).^2-(9*y-7).^2);
end

function C = piston(W,S,V0,k,P0,Ta,T0)
% This function is found in [Ben-Ari and Steinberg (2007)]
A = P0.*S + 19.62*W - k.*V0./S;
PVTT = P0.*V0.*Ta./T0;
V = S./(2*k).*(sqrt(A.^2 + 4*k.*PVTT) - A);
C = 2*pi*sqrt(W./(k+(S./V).^2.*PVTT));
end
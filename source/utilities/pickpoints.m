function [x,spacestr] = pickpoints(a,b,N,spaceopt,ep,alpha)
% function [x,spacestr] = pickpoints(a,b,N,spaceopt,ep,alpha)
% This function creates points in 1D according to a spacing you choose
%
% Inputs - a : the lower bound for the points
%          b : the upper bound for the points
%          N : number of points requested
%          spaceopt : how you want the points spaced
%             even : <default> evenly spaced points
%             cheb : Chebyshev nodes
%             inner : Gaussian cluster near the center
%             halton : described by the Halton sequence
%          ep : the Gaussian shape parameter
%             Note that ep is only needed for 'inner' and 'roots' points
%          alpha : the GaussQR global scale parameter
%             Note that alpha is only need for 'roots' points
%
% Outputs - x : a column vector of points spaced as you requested
%           spacestr : a string describing the spacing of the points
if nargin<6
    alpha = 1;
    if nargin<5
        ep=0;
        if nargin<4
            spaceopt='even';
        end
    end
end

switch lower(spaceopt)
    case 'even'
        x = linspace(a,b,N)';
        spacestr=' evenly spaced';
    case 'inner'
        sigma = 1/(2*sqrt(ep^2+.5*sqrt(4*ep^4+1)));
        x = icdf('Normal',linspace(1/(2*N-1),(2*N-1)/(2*N),N),0,sigma)';
        xmax = max(x);xmin = min(x);xscale = (b-a)/(xmax-xmin);xpt = a-xscale*xmin;
        x = xpt + xscale*x;
        spacestr=' central cluster';
    case 'cheb'
        x = -.5*(b-a)*cos(pi*(0:N-1)/(N-1))' + .5*(b+a);
        spacestr=' Chebyshev points';
    case 'halton'
        x = (b-a)*haltonseq(N,1) + a;
        spacestr=' halton points';
    case 'rand'
        x = (b-a)*rand(N,1) + a;
        spacestr=' random points';
    case 'roots'
        x = gqr_roots(N+1,ep,alpha,1);
        xmax = max(x);
        x = (b-a)/(2*xmax)*x + (b+a)/2;
        spacestr=' eigenfunction points';
end
function [x,spacestr] = pickpoints(a,b,N,spaceopt,ep)
% Note that ep is only needed for 'inner' points
if nargin<5
    ep=0;
    if nargin<4
        spaceopt='even';
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
        x = -.5*(b-a)*cos(pi*(0:N-1)/(N-1))'+.5*(b+a);
        spacestr=' Chebyshev points';
    case 'halton'
        x = (b-a)*haltonseq(N,1)+a;
        spacestr=' halton points';
end
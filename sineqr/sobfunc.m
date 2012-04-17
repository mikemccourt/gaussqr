function [s,M] = sobfunc(x,z,L,sigma,beta,deriv,Mfix)
% function [s,M] = sobfunc(x,z,L,sigma,beta,deriv,Mfix)
% This function evaluates the kernel associated with:
%   something that I need to fill in later
% We can only evaluate it with this summation
% This function will be 0 at x=0 and x=L
%
% When this function is called with beta=1 and Mfix=0, it will use the
% closed form we have available.  If passed Mfix, it will override the
% closed form with the series approximation.
%
% Note that if you pass an x value too long, this probably will crash
% because it will run out of memory.  If that happens, please call it in
% chunks, which will allow it to evaluate successfully
%
% Inputs : x - column vector of evaluation points
%          z - column vector of kernel centers
%            * IF z is a single value, this will use that as all centers
%          L - size of domain of BVP
%          sigma - peakedness parameter of kernel
%          beta - <optional, default = 1> order of the differential operator
%          deriv - <optional, default = 0> order of derivative
%          Mfix - <optional> demanded length of the summation
% Outputs : s - kernel function evaluation
%           M - <optional> length of summation used
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
sumtol = GAUSSQR_PARAMETERS.SUMMATION_TOLERANCE;

% This allows for passing just one center if it is the same for all the
% functions that we want to evaluate
if length(z)==1
    z = z*ones(size(x));
end
if nargin<7
    Mfix = 0;
    if nargin<6
        deriv = 0;
        if nargin<5
            beta = 1;
        end
    end
end

[rx,cx] = size(x);
[rz,cz] = size(z);

if rx~=rz
    error('dimension mismatch, size(x,1)=%d, size(z,1)=%d',rx,rz)
elseif cx~=1 | cz~=1
    error('only column vectors for x and z allowed, size(x,2)=%d, size(z,2)=%d',cx,cz)
elseif L<=0 | length(L)~=1 | imag(L)~=0
    error('unacceptable length value L=%g',L)
elseif sigma<=0 | length(sigma)~=1 | imag(sigma)~=0
    error('unacceptable peakedness parameter sigma=%g',sigma)
end

if beta<1 | length(beta)~=1 | imag(beta)~=0 | floor(beta)~=beta
    warning('unacceptable order beta=%g, reset to 1',beta)
end
if deriv<0 | deriv>4 | length(deriv)~=1 | imag(deriv)~=0 | floor(deriv)~=deriv
    warning('unacceptable derivative deriv=%g, reset to 0',deriv)
    deriv = 0;
end
if Mfix<0 | length(Mfix)~=1 | imag(Mfix)~=0 | floor(Mfix)~=Mfix
    warning('unacceptable length Mfix=%g, choosing internally',Mfix)
    Mfix = 0;
end

if deriv>0 & beta==1
    error('First order kernel is not differentiable, deriv=%d, beta=%d',deriv,beta)
elseif deriv>2 & beta==2
    error('Second order kernel can have at most 2 derivatives, deriv=%d, beta=%d',deriv,beta)
end

% These are the eigenfunctions of the series
% Each line below is a derivative
switch deriv
    case 0
        sinfunc = @(n,L,x) sqrt(2/L)*sin(pi*x*n/L);
    case 1
        sinfunc = @(n,L,x) (sqrt(2/L)*ones(size(x))*(pi*n/L)).*cos(pi*x*n/L);
    case 2
        sinfunc = @(n,L,x) (-sqrt(2/L)*ones(size(x))*(pi*n/L).^2).*sin(pi*x*n/L);
    case 3
        sinfunc = @(n,L,x) (-sqrt(2/L)*ones(size(x))*(pi*n/L).^3).*cos(pi*x*n/L);
    case 4
        sinfunc = @(n,L,x) (sqrt(2/L)*ones(size(x))*(pi*n/L).^4).*sin(pi*x*n/L);
    otherwise
        error('Unacceptable derivative called, deriv=%d',deriv)
end
% These are the eigenvalues of the series
lamfunc = @(n,L,sigma,beta) ((pi*n/L).^2+sigma^2).^(-beta);

if Mfix==0
    M = floor(1/pi*sqrt(sumtol^(-1/beta)*(pi^2+(sigma*L)^2)-(sigma*L)^2));
else
    M = Mfix;
end

% For beta==1 we have the closed form of the function
% Otherwise we need to compute the series
if beta==1 && Mfix==0
    s = sinh(sigma*min(x,z)).*sinh(sigma*(L-max(x,z)))./(sigma*sinh(L*sigma));
else
    % Should come up with a way to automate this so that it will check for
    % allocation issues and work around it.  Maybe it could be optimized to
    % account for cache misses?
    blocksize = 1000;
    s = zeros(size(x));

    Mind = M;
    for k=1:floor(M/blocksize)
        n = Mind-blocksize+1:Mind;

        Xmat = sinfunc(n,L,x);
        Zmat = sinfunc(n,L,z);
        Xmat = Xmat.*Zmat;
        Zmat = repmat(lamfunc(n,L,sigma,beta),rx,1);
        Xmat = Xmat.*Zmat;

        s = s + sum(Xmat,2);

        Mind = Mind - blocksize;
    end
    if Mind>0
        n = 1:Mind;

        Xmat = sinfunc(n,L,x);
        Zmat = sinfunc(n,L,z);
        Xmat = Xmat.*Zmat;
        Zmat = repmat(lamfunc(n,L,sigma,beta),rx,1);
        Xmat = Xmat.*Zmat;

        s = s + sum(Xmat,2);
    end
end



end
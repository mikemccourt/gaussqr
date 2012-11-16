function [L,ep,beta,Mmax] = mqr_solveprep(x,L,ep,beta)
% function [L,ep,beta,Mmax] = mqr_solveprep(x,L,ep,beta)
%
% This function takes all the stuff needed to do a GaussQR problem and
% preps it before the problem starts.  The reason for this is so that I can
% get some meta-info without calling the whole solve routine.  Also,
% changes which are the same for both QR and QRr will only need changed at
% one place rather than multiple files
%
% INPUTS: x - input data values
%         L - boundary of the domain
%         ep - value of shape parameter
%         beta - <default=1> smoothness parameter
% OUTPUTS: L - boundary of the domain
%          ep - value of shape parameter
%          beta - smoothness parameter
%          Mmax - length of the necessary eigenfunction expansion
%
% NOTE: When this function is called with nothing passed, it will return a
% function handle allowing you to form the eigenvalues.  That function has
% the calling structure
%    lamfunc = @(n,L,ep,beta) ((pi*n/L).^2+ep^2).^(-beta);
% n is a row vector of indices

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
Mextramax = GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC;

if nargin==0
    L = @(n,L,ep,beta) ((pi*n/L).^2+ep^2).^(-beta);
else
    [N,d] = size(x);

    % Checks to make sure that the ep and beta values are acceptable
    if length(ep)>1
        ep = abs(real(ep(1)));
        warning('Multiple epsilon values not allowed; using epsilon=%g',ep)
    end
    if abs(real(ep))~=ep
        ep = abs(real(ep));
        warning('Only real, positive epsilon allowed; using epsilon=%g',ep)
    end
    if length(beta)>1
        beta = abs(real(beta(1)));
        warning('Multiple beta values not allowed; using beta=%g',beta)
    end
    if beta~=floor(abs(real(beta)))
        beta = floor(abs(real(beta)));
        warning('Non-integer beta values not yet allowed; using beta=%d',beta)
    end

    % Checks to make sure that L is an acceptable value
    if x~=abs(real(x))
        error('x must have nonnegative, real values only')
    end
    minx = min(x);
    maxx = max(x);
    if minx<0
        error('Values for MaternQR must be on [0,L], min(x)=%g',minx)
    elseif L==maxx
        warning('Are you SURE you meant to pass x=0 or x=L?')
    elseif L<maxx
        L = maxx;
        warning('Data passed beyond [0,L] domain, resetting L=%g',L)
    end

    % This picks the eigenfunction series length that we want to use for the
    % mqr_solve.  The minimum length is 1.01*N, although eventually I'll allow
    % for smaller M values
    if Mextramax<0
        Mmax = ceil((1-(min(Mextramax,-101)/100))*N);
    else
        Mmax = ceil(max(Mextramax,1.01*N));
    end

end

function [ep,alpha,Marr,lam] = gqr_solveprep(reg,x,ep,alpha,M)
% function [ep,alpha,Marr,lam] = gqr_solveprep(reg,x,ep,alpha,M)
%
% This function takes all the stuff needed to do a GaussQR problem and
% preps it before the problem starts.  The reason for this is so that I can
% get some meta-info without calling the whole solve routine.  Also,
% changes which are the same for both QR and QRr will only need changed at
% one place rather than multiple files
%
% INPUTS: reg - pass a 1 for regression, 0 otherwise
%         x - input data values
%         ep - value of epsilon shape parameter
%         alpha - (optional) value of global scale parameter
%         M - (optional) truncation value suggested by user
%
% Note: for regression, there are only 3 outputs (not lam)
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
Mextramax = GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC;
Mdefault = GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC;
alphaDefault = GAUSSQR_PARAMETERS.ALPHA_DEFAULT;
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;

[N,d] = size(x);

% If the user didn't pass alpha, or passed [], then we need to pick an
% alpha from rbfalphasearch
computealpha = 0;
if nargin==3
    computealpha = 1;
elseif length(alpha)==0
    computealpha = 1;
end
if computealpha==1
    xminBound = min(x);
    xmaxBound = max(x);
    alpha = gqr_alphasearch(ep,xminBound,xmaxBound);
end

% Checks to make sure that the ep and alpha values are acceptable
if length(ep)>1
    ep = abs(real(ep(1)));
    warning(sprintf('Multiple epsilon values not allowed; using epsilon=%g',ep))
end
if length(alpha)>1
    alpha = abs(real(alpha(1)));
    warning(sprintf('Multiple alpha values not allowed; using alpha=%g',alpha))
end
if abs(real(ep))~=ep
    ep = abs(real(ep));
    warning(sprintf('Only real, positive epsilon allowed; using epsilon=%g',ep))
end
if abs(real(alpha))~=alpha
    alpha = abs(real(alpha));
    warning(sprintf('Only real, positive alpha allowed; using alpha=%g',alpha))
end
if ep==0 || alpha==0
    error(sprintf('Parameters cannot be zero: epsilon=%g, alpha=%g',ep,alpha))
end

% This switch changes what we're computing based on if the
% user wants to do interpolation or regression
if reg
    % If we call this function but don't actually want to compute Marr we
    % can skip the rest of this
    if nargout==3
        if Mdefault<=1 % Only a percentage of N is considered
            Mdefault = floor(Mdefault*N);
        elseif Mdefault>N
            Mdefault = N;
        end

        % All this figures out what an acceptable Marr is
        % For now I'm doing Marr+1 to solve for an off-by-one issues in the paper
        % I'll work on fixing this eventually
        if not(exist('M'))
            M = zeros(d,1);
            Mlim = Mdefault;
        else
            % Need to check to make sure passed M is reasonable
            if any(M>N)
                warning('M value %d>%d unacceptable, reset to %d',M,N,Mdefault)
                Mlim = Mdefault;
                M = zeros(d,1);
            elseif any(M<0)
                warning('Negative M value passed, setting max Marr length to %d',Mdefault)
                Mlim = Mdefault;
                M = zeros(d,1);
            end
            [Md,Mc] = size(M);
            if Md~=d
                if Md==1 % The user only passed a max length for Marr
                    Mlim = M;
                    M = zeros(d,1);
                else
                    error('M of length %d passed, but data has dimension %d',Md,d)
                end
            elseif Mc~=1
                error('Multiple columns in M detected, but only 1 is allowed')
            else % This is the everything was passed correctly case
                Mlim = N;
            end
        end

        Marr = gqr_formMarr(M,[],Mlim);
    end
else
    nu = (2*ep/alpha)^2;
    lam = nu/(2+nu+2*sqrt(1+nu));
    if Mextramax<0
        Mextramax = (1-Mextramax/100)*N;
    end
    MarrN = gqr_formMarr(zeros(d,1),[],N);
    Mlim = ceil(size(MarrN,2)+log(eps)/log(lam));
    Mlim = ceil(N+log(eps)/log(lam));
    if Mextramax==0
        Mextramax = inf; % Allow the array to go as long as it wants
    end

    % This needs to get better
    % Specifically it needs to handle people passing weird stuff
    if not(exist('M'))
        M = zeros(d,1);
    else
        [Mr Mc] = size(M);
        if Mr~=d
            error('Incorrect M size passed, size(M)=%dx%d d=%d',Mr,Mc,d)
        elseif M<N
            error('gqr_solve requires M>N, but M=%g, N=%d',M,N)
        elseif ceil(M)~=M
            warning('Noninteger M passed as %g, reset to %d',M,ceil(M))
            M = ceil(M);
        end
    end

    Marr = gqr_formMarr(M,Mlim,Mextramax);
end
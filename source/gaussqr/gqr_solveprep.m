function [ep,alpha,Marr,Rbar] = gqr_solveprep(reg,x,ep,alpha,M)
% This function takes all the stuff needed to do a GaussQR problem and
% preps it before the problem starts.  The reason for this is so that I can
% get some meta-info without calling the whole solve routine.  Also,
% changes which are the same for both QR and QRr will only need changed at
% one place rather than multiple files
%
% function GQR = gqr_solveprep(reg,x,ep)
% Inputs : reg - pass 1 for regression, 0 interpolation
%          x - input data values
%          ep - value of epsilon shape parameter
%          alpha - (optional) value of global scale parameter
%          M - (optional) truncation value suggested by user
% Outputs : GQR - the GaussQR object, with the needed members
%   NOTE: reg=-1 allows for creating a GQR shell and testing the input x,
%   ep, alpha values for validity.  This is nice if you want to compute the
%   eigenvalues of certain ep and alpha values with GQR.eig()
%
% function GQR = gqr_solveprep(reg,x,ep,alpha)
% Inputs : alpha - (optional) value of global scale parameter
%
% function GQR = gqr_solveprep(reg,x,ep,alpha,M)
% Inputs : M - (optional) truncation value suggested by user
%     If you want to pass an M, but not an alpha, use
%          GQR = gqr_solveprep(reg,x,ep,[],M)
%
% function [ep,alpha,Marr,Rbar] = gqr_solveprep(0,...)
% Outputs : ep - the acceptable shape parameter
%           alpha - the acceptable scale parameter
%           Marr - the GQR index list for interpolation
%           Rbar - the matrix such that psi = gqr_phi*[I;Rbar]
%
% function [ep,alpha,Marr] = gqr_solveprep(1,...)
% Outputs : ep - the acceptable shape parameter
%           alpha - the acceptable scale parameter
%           Marr - the GQR index list for regression
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
Mextramax = GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC;
Mdefault = GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC;
alphaDefault = GAUSSQR_PARAMETERS.ALPHA_DEFAULT;
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;
storephi = GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION;
asympttol = GAUSSQR_PARAMETERS.RBFPHI_EXP_TOL;

computealpha = 0;
switch nargin
    case {0,1,2}
        error('nargin=%d is too few input arguments',nargin)
    case 3
        computealpha = 1;
    case {4,5}
        if length(alpha)==0;
            computealpha = 1;
        end
    otherwise
        error('nargin=%d is too large',nargin)
end

returnGQR = 0;
switch nargout
    case 1
        returnGQR = 1;
    case 3
        if reg==0
            error('For reg=%g (interpolation), 3 outputs is unacceptable',reg)
        end
    case 4
        if reg==1
            error('For reg=%g (regression), 4 outputs is unacceptable',reg)
        end
    otherwise
        error('nargout=%d unacceptable',nargout)
end

[N,d] = size(x);

% If the user didn't pass alpha, or passed [], then we need to pick an
% alpha from gqr_alphasearch
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

% Evaluate auxiliary parameters
beta = (1+(2*ep/alpha)^2)^(1/4);
if beta-1<asympttol % This triggers an asymptotic expansion
    delta2 = ep^2-ep^4/alpha^2+2*ep^6/alpha^4;
else
    delta2 = 1/2*alpha^2*(beta^2-1);
end

% Set up GQR object for evaluating eigenfunctions if necessary
GQR.reg = reg;
GQR.ep = ep;
GQR.alpha = alpha;
GQR.beta = beta;
GQR.delta2 = delta2;
GQR.eig = @(n) alpha/sqrt(alpha^2+ep^2+delta2)*(ep^2/(alpha^2+ep^2+delta2)).^(sum(n,1)-d);
GQR.warnid = '';
GQR.warnmsg = '';
if storephi
    GQR.stored_x = [];
end

% This switch changes what we're computing based on if the
% user wants to do interpolation or regression
switch reg
    case 1
        if Mdefault<=1 % Only a percentage of N is considered
            Mdefault = floor(Mdefault*N);
        elseif Mdefault>N
            Mdefault = N;
        end

        % All this figures out what an acceptable Marr is
        if not(exist('M','var'))
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
        GQR.Marr = Marr;
    case {0,-1} % Create the GQR shell or the whole interpolation object
        % First create Marr
        % We need the eigenvalues for this
        % Only the relative eigenvalue ratio is needed
        % Thus the constant alpha/sqrt(ep^2+alpha^2+delta2) is omitted
        lam = ep^2/(ep^2+alpha^2+delta2);
        
        if Mextramax<0
            Mextramax = (1-Mextramax/100)*N;
        elseif Mextramax>0
            Mextramax = N + Mextramax;
        else
            Mextramax = inf; % Allow the array to go as long as it wants
        end

        MarrN = gqr_formMarr(zeros(d,1),[],N);
        Mlim = ceil(N+log(eps)/log(lam));

        % This needs to get better
        % Specifically it needs to handle people passing weird stuff
        if not(exist('M','var'))
            M = zeros(d,1);
        else
            [Mr,Mc] = size(M);
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
        GQR.Marr = Marr;
        
        if reg==0 % Stop now if only the shell was requested
            % Now we need to create the Rbar matrix
            phiMat = gqr_phi(GQR,x);
            [Q,R] = qr(phiMat);
            R1 = R(:,1:N);
            R2 = R(:,N+1:end);
            
            lastwarn('')
            warning off MATLAB:divideByZero
            iRdiag = diag(1./diag(R1));
            [warnmsg,msgid] = lastwarn;
            if strcmp(msgid,'MATLAB:divideByZero')
                GQR.warnid = 'GAUSSQR:zeroQRDiagonal';
                GQR.warnmsg = 'At least one value on the R diagonal was exactly 0';
            end
            warning on MATLAB:divideByZero
            
            R1s = iRdiag*R1;
            opts.UT = true;
            
            lastwarn('')
            warning off MATLAB:singularMatrix
            Rhat = linsolve(R1s,iRdiag*R2,opts);
            warning on MATLAB:singularMatrix
            
            [warnmsg,msgid] = lastwarn;
            warnid = 'GAUSSQR:singularR1invR2';
            warnmsg = 'Computing inv(R1)R2 ... R1 singular to working precision';
            
            if ~strcmp(GQR.warnid,'GAUSSQR:zeroQRDiagonal') && strcmp(msgid,'MATLAB:singularMatrix')
                GQR.warnid = warnid;
                GQR.warnmsg = warnmsg;
            end
            
            % Here we apply the eigenvalue matrices
            % Note that the -d term in the power goes away because it
            % appears in both Lambda_2 and Lambda_1^{-1}
            Ml = size(Marr,2);
            D = lam.^(repmat(sum(Marr(:,N+1:end),1)',1,N)-repmat(sum(Marr(:,1:N),1),Ml-N,1));
            Rbar = D.*Rhat';
            
            GQR.Rbar = Rbar;
            if storephi && returnGQR
                GQR.stored_x = x;
                GQR.stored_deriv = 0;
                GQR.stored_phi1 = phiMat(:,1:N);
                GQR.stored_phi2 = phiMat(:,N+1:end);
            end
        end
    otherwise
        error('reg=%g is unacceptable, 1 for regression, 0 for interpolation, -1 for eigenvalue setup',reg)
end
        
if returnGQR
    ep = GQR;
end

if alertuser && ~strcmp(GQR.warnid,'')
    warning(GQR.warnid,GQR.warnmsg)
end

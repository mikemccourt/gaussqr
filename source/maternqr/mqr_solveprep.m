function [L,ep,beta,M,Rbar] = mqr_solveprep(x,L,ep,beta)
% function [L,ep,beta,M,Rbar] = mqr_solveprep(x,L,ep,beta)
%
% This function takes all the stuff needed to do a MaternQR problem and
% preps it before the problem starts.  The reason for this is so that I can
% get some meta-info without calling the whole solve routine.
%
% function lamfunc = mqr_solveprep()
% Outputs : lamfunc -  a function handle to evaluat the eigenvalues.
%                      lamfunc has the calling structure
%         lamfunc = @(n,L,ep,beta) ((pi*n/L).^2+ep^2).^(-beta);
%                      n is a row vector of indices
%
% function MQR = mqr_solveprep(x,L,ep,beta)
% Inputs : x - input data values
%          L - boundary of the domain
%          ep - value of shape parameter
%          beta - smoothness parameter (integer >=1)
% Outputs : MQR - the MathernQR object, with the needed members
%
% function [L,ep,beta,Mmax,Rbar] = mqr_solveprep(x,L,ep,beta)
% Outputs: L - boundary of the domain
%          ep - value of shape parameter
%          beta - smoothness parameter
%          M - length of the necessary eigenfunction expansion
%          Rbar - the stable basis correction coefficients

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
Mextramax = GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC;
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;
storephi = GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION;

% Define the function which evaluates eigenvalues
lamfunc = @(n,L,ep,beta) ((pi*n/L).^2+ep^2).^(-beta);

if nargin==0
    L = lamfunc;
elseif nargin~=4
    error('nargin=%d is unacceptable',nargin)
else
    returnMQR = 0;
    if nargout==1
        returnMQR = 1;
    end
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
        M = ceil((1-(min(Mextramax,-101)/100))*N);
    else
        M = ceil(max(Mextramax,1.01*N));
    end
    
    MQR.ep = ep;
    MQR.L = L;
    MQR.beta = beta;
    MQR.N     = N;
    MQR.warnid = '';
    MQR.warnmsg = '';
    if storephi
        MQR.stored_x = [];
    end
    
    % Define the eigenvalues
    lamfunc = mqr_solveprep();
    
    % Stuff needed for the system solve
    I = eye(N);
    opts.UT = true;
    
    % Form the Marr for this system, just 1:Mmax
    Marr = 1:M;
    phiMat = mqr_phi(Marr,x,L);
    
    % Compute the QR decomposition of the short, fat matrix
    [Q,R] = qr(phiMat);
    R1 = R(:,1:N);
    R2 = R(:,N+1:end);
    
    lastwarn('')
    
    warning off MATLAB:singularMatrix
    warning off MATLAB:nearlySingularMatrix
    
    % Apply inv(R1)*R2
    [Rhat,rcond] = linsolve(R1,R2,opts);
    [warnmsg,msgid] = lastwarn;
    if strcmp(msgid,'MATLAB:singularMatrix')
        MQR.warnid = 'MATERNQR:singularR1invR2';
        MQR.warnmsg = 'Computing inv(R1)R2 ... R1 singular to working precision';
    end
    if strcmp(msgid,'MATLAB:nearlySingularMatrix')
        MQR.warnid = 'MATERNQR:nearlySingularR1invR2';
        MQR.warnmsg = 'Computing inv(R1)R2 ... R1 nearly singular, rcond = %g',rcond;
    end
    
    warning off MATLAB:singularMatrix
    warning off MATLAB:nearlySingularMatrix
    
    % Compute the eigenvalue matrix
    lambda = lamfunc(Marr,L,ep,beta);
    D1 = repmat(lambda(1:N),M-N,1);
    D2 = repmat(lambda(N+1:end)',1,N);
    
    % Form the Rbar matrix
    Rbar = D2.*Rhat'./D1;
    
    MQR.Marr  = Marr;
    MQR.Rbar  = Rbar;
    
    if storephi
        MQR.stored_x = x;
        MQR.stored_deriv = 0;
        MQR.stored_phi1 = phiMat(:,1:N);
        MQR.stored_phi2 = phiMat(:,N+1:end);
    end
    
    if alertuser && ~strcmp(MQR.warnid,'')
        warning(MQR.warnid,MQR.warnmsg)
    end
    
    if returnMQR
        L = MQR; % Store it in the right object to return
    end
    
end

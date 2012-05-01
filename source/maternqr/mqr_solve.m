function MQR = mqr_solve(x,y,L,ep,beta)
% function MQR = mqr_solve(x,y,L,ep,beta)
% This function finds the interpolant to scattered data in 1D using the
% eigenfunction expansion of the Compact Matern kernel
% The problem is defined on the domain [0,L]
%
% Inputs: x - column vector of data points between [0,L]
%         y - column vector of data values evaluated at x
%         L - boundary of the domain, must be positive
%         ep - shape parameter, must be positive
%         beta - <default=1> smoothness order of spline
% Outputs: MQR - maternQR object with relevant data stored

% Import global parameters from rbfsetup
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;

if nargin<4
    error('Insufficient inputs')
end
MQR.warnid = '';
MQR.warnmsg = '';

[N,d] = size(x);
if d~=1
    error('Only 1D problems are acceptable for MaternQR, size(x,2)=%d',d)
elseif sum(N~=size(y,1))
    error('Different number of input and output points (x and y)')
elseif size(y,2)~=1
    error('You can only pass a 1D output vector y')
end

if nargin==4
    [L,ep,beta,M] = mqr_solveprep(x,L,ep);
else
    [L,ep,beta,M] = mqr_solveprep(x,L,ep,beta);
end

% Define the eigenvalues
lamfunc = mqr_solveprep();

% Stuff needed for the system solve
I = eye(N);
opts.UT = true;

% Form the Marr for this system, just 1:Mmax
n = 1:M;
S = mqr_phi(n,x,L);
% S = sinfunc(n,L,x);

% Compute the QR decomposition of the short, fat matrix
[Q,R] = qr(S);
R1 = R(:,1:N);
R2 = R(:,N+1:end);

lastwarn('')

warning off MATLAB:singularMatrix
warning off MATLAB:nearlySingularMatrix

% Apply inv(R1)*R2
[Rhat,rcond] = linsolve(R1,R2,opts);
[warnmsg,msgid] = lastwarn;
if strcmp(msgid,'MATLAB:singularMatrix')
    MQR.warnid = 'GAUSSQR:singularR1invR2';
    MQR.warnmsg = 'Computing inv(R1)R2 ... R1 singular to working precision';
end
if strcmp(msgid,'MATLAB:nearlySingularMatrix')
    MQR.warnid = 'GAUSSQR:nearlySingularR1invR2';
    MQR.warnmsg = 'Computing inv(R1)R2 ... R1 nearly singular, rcond = %g',rcond;
end

% Compute the eigenvalue matrix
lambda = lamfunc(n,L,ep,beta);
D1 = repmat(lambda(1:N),M-N,1);
D2 = repmat(lambda(N+1:end)',1,N);

% Form the Rbar matrix
Rbar = D2.*Rhat'./D1;

lastwarn('')

% Solve for the interpolant coefficients
[coef,rcond] = linsolve(S*[I;Rbar],y);
[warnmsg,msgid] = lastwarn;
if strcmp(msgid,'MATLAB:singularMatrix')
    MQR.warnid = 'GAUSSQR:singularPsiSolve';
    MQR.warnmsg = 'Computing Psi\y singular to working precision';
end
if strcmp(msgid,'MATLAB:nearlySingularMatrix')
    MQR.warnid = 'GAUSSQR:nearlySingularPsiSolve';
    MQR.warnmsg = 'Computing Psi\y nearly singular, rcond = %g',rcond;
end

warning on MATLAB:singularMatrix
warning on MATLAB:nearlySingularMatrix

% Populate the MQR object with the necessary values
MQR.ep    = ep;
MQR.L     = L;
MQR.beta  = beta;
MQR.N     = N;
MQR.Mmax  = M;
MQR.Rbar  = Rbar;
MQR.coef  = coef;

if alertuser && ~strcmp(MQR.warnid,'')
    warning(MQR.warnid,MQR.warnmsg)
end
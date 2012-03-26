function eigroots = rbfroots(M,ep,alpha,allroots)
% function eigroots = rbfroots(M,ep,alpha,allroots)
% This function returns the roots of the eigenfunction defined by
%    @(x) rbfphi(M,x,ep,alpha)
% It does so by forming a symmetric tridiagonal matrix and using its
% eigenvalues as the roots
%
% Inputs : M - order of the eigenfunction (M>=2)
%          ep - Gaussian shape parameter
%          alpha - GaussQR scale parameter
%          allroots - (optional) how many roots you want
%                     generally you only want the largest root,
%                     but pass allroots=1 to get all roots
%
% Outputs: eigroots - the roots you requested, in increasing order

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
splimit = GAUSSQR_PARAMETERS.ORTH_ROOT_SPARSE_LIMIT;

if M<2
    error('First eigenfunction has no roots')
end

if nargin==3
    allroots = 0;
end

if M<splimit | allroots==1
    J = diag(sqrt(1:M-2),1)+diag(sqrt(1:M-2),-1);
    eigvals = sort(eig(J));J,eig(J)
    if allroots==1
        eigroots = eigvals;
    else
        eigroots = eigvals(end);
    end
else
    sM = spdiags(sqrt(1:M-2)',0,M-1,M-1);
    Jf = @(x) [zeros(1,size(x,2));sM*x(1:end-1,:)] + [sM*x(2:end,:);zeros(1,size(x,2))];
    opts.display = 0;
    eigroots = eigs(Jf,M-1,1,'lm',opts);
end

% Finish by scaling the roots according to the domain
eigroots = eigroots/(sqrt(2)*alpha*(1+(2*ep/alpha)^2)^.25);
%function errr = errorevals(Phi,M)
%Evaluates the error for the eigenvalue
%
%Inputs : phi - The object get form HSeigsolve
%
%Output : errr - a vector of eigenvalues' relative error
function errr = errorevals(Phi)
    % compute eigenvalues' relative error
    if (isfield(Phi,'epsilon'))
        eps = Phi.epsilon;
    else
        eps = 0;
    end
    realr = @(n) 1./((n.^2)*pi^2+eps^2);
    D = Phi.eigvals;
    s = size(D);
    n = [1: s(1)]';
    Rer = realr(n);
    errr = abs( D - Rer)./Rer;
  
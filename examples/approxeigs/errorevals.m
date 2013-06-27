%function errr = errorevals(Phi,M)
%Evaluates the error for the eigenvalue
%
%Inputs : phi - The object get form HSeigsolve
%
%Output : errr - a vector of eigenvalues' relative error
function errr = errorevals(Phi)
    % compute eigenvalues' relative error 
    realr = @(n) 1./((n.^2)*pi^2);
    D = Phi.eigvals;
    C = Phi.coefs;
    s = size(D);
    n = [1: s(1)]';
    Rer = realr(n);
    errr = abs( D - Rer)./Rer;
  
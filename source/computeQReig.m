function [invU,Svec,Q] = computeQReig(M,x,ep,alpha)
% function [invU,Svec,Q] = computeQReig(M,x,ep,alpha)
%
% This computes the decomposition
%    P = Q * S * U
% where: P - rbfphi(Mx,ep,alpha)
%        Q - from P = QR, has orthonormal columns
%        S - diagonal matrix, such that U'*S^2*U = P'*P
%        U - from the LDL decomposition of U'*S^2*U = P'*P
%
% Inputs: M - length of the eigenfunction summation
%             note that M<size(x,1)
%         x - the collocation points
%         ep - Gaussian shape parameter
%         alpha - GaussQR global scale parameter
%
% Outputs: invU - inv(U) where U is defined above
%          Svec - diag(S) where S is defined above
%          Q - the orthonormal matrix from the QR decomposition
%
% Note that this may have less stability than the traditional QR
% factorization via Householder or whatever.  I still have some more
% stability analysis to do, but you should be aware that for larger M it
% could be an issue.
%
% NOTE: TM should have a 1/(beta*alpha) term in it, but that is applied
% analytically throughout because TM always appears with a beta*alpha.  It
% seemed unnecessary to apply it and then immediately undo it.
%
% Possible optimization may exist for computing b.  In the Hermite paper,
% it was computed without an inner product, but I haven't been able to
% determine how that was done, so I haven't been able to make it happen
% here.  Will look at that soon.

N = size(x,1);
invU = zeros(M);
Svec = zeros(M,1);
Q = zeros(N,M);

Dx = (1+(2*ep/alpha)^2)^.25*alpha*x;

phi = rbfphi(1:M,x,ep,alpha);

tmp = phi'*phi(:,[1,M-1:M]);
d = tmp(:,2)*sqrt(M-1)/sqrt(2) - ApplyTM(tmp(:,3));

% Start with the first columnns
u0 = eye(M,1);
v0 = tmp(:,1);

Svec(1) = sqrt(u0'*v0);
invU(:,1) = u0;

q0 = phi(:,1)/Svec(1);
Q(:,1) = q0;

% Get the second columns
tmp = ApplyTM(u0);
c = -tmp'*v0/Svec(1)^2;
u1 = sqrt(2)*(tmp+c*u0);
v1 = sqrt(2)*(ApplyTM(v0) + c*v0);
v1(M) = v1(M) - sqrt(2)*(d'*u0);

Svec(2) = sqrt(u1'*v1);
invU(:,2) = u1;

q1 = (sqrt(2)*Svec(1)/Svec(2))*(Dx.*q0 + c*q0);
Q(:,2) = q1;

% Iterate through the remaining columns
for k = 2:M-1
    tmp = ApplyTM(u1);
    c = -(v1'*tmp)/Svec(k)^2;
    b = .5*(k-1)*Svec(k)^2/Svec(k-1)^2;
    u2 = sqrt(2/k)*(tmp+c*u1) - (2*b/sqrt(k^2-k))*u0;
    v2 = sqrt(2/k)*(ApplyTM(v1)+c*v1) - (2*b/sqrt(k^2-k))*v0;
    v2(M) = v2(M) - sqrt(2/k)*(d'*u1);
    
    Svec(k+1) = sqrt(u2'*v2);
    invU(:,k+1) = u2;
    
    q2 = (sqrt(2/k)*Svec(k)/Svec(k+1))*(Dx.*q1+c*q1) - (2*b/sqrt(k^2-k)*Svec(k-1)/Svec(k+1))*q0;
    Q(:,k+1) = q2;

    u0 = u1; u1 = u2;
    v0 = v1; v1 = v2;
    q0 = q1; q1 = q2;
end
end

% This is a private function that shouldn't be accessed outside this file
% It applies the matrix
%      TM = (diag(sqrt(1:M-1),1) + diag(sqrt(1:M-1),-1))/sqrt(2)
% to a vector without forming that matrix
% M is the length of the vector
%
% NOTE: this function avoids being passed ba by applying it analytically
% throughout above.  Therefore, this is NOT the same TM matrix that appears
% in the paper.
%
% Also, it only should be passed column vectors, not multiple columns.
% This shouldn't be a problem, but I'm leaving a note in case I go senile.
function TMx = ApplyTM(x)
persistent M Mvec k;
if isempty(Mvec)
    M = size(x,1);
    k = 1;
    Mvec = sqrt(1:M-1)'/sqrt(2);
end
TMx = ([Mvec.*x(2:end);0] + [0;Mvec.*x(1:end-1)]);
end
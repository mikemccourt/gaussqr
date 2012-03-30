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

N = size(x,1);
invU = zeros(M);
Svec = zeros(M,1);
Qvec = zeros(M,1);
Q = zeros(N,M);

beta = (1+(2*ep/alpha)^2)^.25;
ba = beta*alpha;
delta2 = .5*alpha^2*(beta^2-1);

phi = rbfphi(1:M,x,ep,alpha);
TM = 1/(sqrt(2)*ba)*(diag(sqrt(1:M-1),-1)+diag(sqrt(1:M-1),1));

tmp = phi'*phi(:,[1,M-1:M]);
d = tmp(:,2)*sqrt(M-1)/(sqrt(2)*ba) - ApplyTM(tmp(:,3),ba);

% Start with the first columnns
u0 = eye(M,1);
v0 = tmp(:,1);

Svec(1) = sqrt(u0'*v0);
invU(:,1) = u0;

q0 = exp(-delta2*x.^2);
q0 = q0/norm(q0);
Q(:,1) = q0;

% Get the second columns

tmp = ApplyTM(u0,ba);
c = -ba*tmp'*v0/Svec(1)^2;
u1 = sqrt(2)*(ba*tmp+c*u0);
v1 = sqrt(2)*(ba*ApplyTM(v0,ba) + c*v0);
v1(M) = v1(M) - sqrt(2)*ba*(d'*u0);

Svec(2) = sqrt(u1'*v1);
invU(:,2) = u1;

q1 = sqrt(2)*(ba*x.*q0 + c*q0)*Svec(1)/Svec(2);
Q(:,2) = q1;

for k = 2:M-1

    tmp = ApplyTM(u1,ba);
    c = -ba/Svec(k)^2*(v1'*tmp);
    b = sqrt((k-1)/2)*ba/Svec(k-1)^2*(v0'*tmp);
    u2 = sqrt(2/k)*(ba*tmp+c*u1) - 2*b/sqrt(k^2-k)*u0;
    v2 = sqrt(2/k)*(ba*ApplyTM(v1,ba)+c*v1) - 2*b/sqrt(k^2-k)*v0;
    v2(M) = v2(M) - sqrt(2/k)*ba*(d'*u1);
    
    Svec(k+1) = sqrt(u2'*v2);
    invU(:,k+1) = u2;
    
    q2 = sqrt(2/k)*(ba*x.*q1+c*q1)*(Svec(k)/Svec(k+1)) - (2*b/sqrt(k^2-k))*(Svec(k-1)/Svec(k+1))*q0;
    Q(:,k+1) = q2;

    u0 = u1; u1 = u2;
    v0 = v1; v1 = v2;
    q0 = q1; q1 = q2;
end
end

% This is a private function that shouldn't be accessed outside this file
% It applies the matrix
%      TM = (diag(sqrt(1:M-1),1) + diag(sqrt(1:M-1),-1))/(sqrt(2)*beta*alpha)
% to a vector without forming that matrix
% M is the length of the vector
%
% NOTE: this function could avoid being passed ba, but I need to make some
% changes throughout to make that viable
%
% Also, it only should be passed column vectors, not multiple columns yet
function TMx = ApplyTM(x,ba)
persistent M Mvec k;
if isempty(Mvec)
    M = size(x,1);
    k = 1;
    Mvec = sqrt(1:M-1)'/(sqrt(2)*ba);
end
TMx = ([Mvec.*x(2:end);0] + [0;Mvec.*x(1:end-1)]);
end
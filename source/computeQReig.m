function [invU,Svec,Q] = computeQReig(M,x,ep,alpha,rhs)
% function [invU,Svec,Q] = computeQReig(M,x,ep,alpha,rhs)
%
% This computes the decomposition
%    P = Q * S * U
% where: P - rbfphi(Mx,ep,alpha)
%        Q - from P = QR, has orthonormal columns
%        S - diagonal matrix, such that U'*S^2*U = P'*P
%        U - from the LDL decomposition of U'*S^2*U = P'*P
%
%%%%%
%
% Calling: [invU,Svec,Q] = computeQReig(M,x,ep,alpha)
%          This computes the full factorization and returns it
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
%%%%%
%
% Calling: sol = computeQReig(M,x,ep,alpha,rhs)
%          This computes only the solution to the problem
%                  Phi*sol = rhs
%          without returning the factorization
%          You can only pass one rhs right now
%
% Inputs: M - length of the eigenfunction summation
%             note that M<size(x,1)
%         x - the collocation points
%         ep - Gaussian shape parameter
%         alpha - GaussQR global scale parameter
%         rhs - right hand side of Phi*sol = rhs
%
% Outputs: sol - column vector such that Phi*sol = rhs
%
% Note that this may have less stability than the traditional QR
% factorization via Householder or whatever.  I still have some more
% stability analysis to do, but you should be aware that for larger M it
% could be an issue.
%
% I should also add a version of this code that takes in a RHS and solves
% the system without actually storing the factorization.

if nargin<4
    error('Insufficient arguments passed')
elseif nargin==4
    [invU,Svec,Q] = computeQRfactorization(M,x,ep,alpha);
else
    invU = computeQRsolve(M,x,ep,alpha,rhs);
    Svec = [];
    Q = [];
end

end

% This computes the actual factorization, and shouldn't be called directly
function [invU,Svec,Q] = computeQRfactorization(M,x,ep,alpha)
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

Svec(1) = sqrt(tmp(1,1));
invU(:,1) = u0;

q0 = phi(:,1)/Svec(1);
Q(:,1) = q0;

% Get the second columns
tmp = ApplyTM(u0);
c = -tmp'*v0/Svec(1)^2;
s2ok = sqrt(2);
u1 = s2ok*(tmp+c*u0);
v1 = s2ok*(ApplyTM(v0) + c*v0);
v1(M) = v1(M) - s2ok*(d'*u0);

Svec(2) = sqrt(u1'*v1);
invU(:,2) = u1;

q1 = (s2ok*Svec(1)/Svec(2))*(Dx.*q0 + c*q0);
Q(:,2) = q1;

% Iterate through the remaining columns
for k = 2:M-1
    tmp = ApplyTM(u1);
    c = -(v1'*tmp)/Svec(k)^2;
    b = .5*(k-1)*(Svec(k)/Svec(k-1))^2;
    s2ok = sqrt(2/k);
    bosk = 2*b/sqrt(k^2-k);
    u2 = s2ok*(tmp+c*u1) - bosk*u0;
    v2 = s2ok*(ApplyTM(v1)+c*v1) - bosk*v0;
    v2(M) = v2(M) - s2ok*(d'*u1);
    
    Svec(k+1) = sqrt(u2'*v2);
    invU(:,k+1) = u2;
    
    q2 = (s2ok*Svec(k)/Svec(k+1))*(Dx.*q1+c*q1) - (bosk*Svec(k-1)/Svec(k+1))*q0;
    Q(:,k+1) = q2;

    u0 = u1; u1 = u2;
    v0 = v1; v1 = v2;
    q0 = q1; q1 = q2;
end
end

% This private function computes the actual solution without storing the
% full factorization.  That allows for less storage
function sol = computeQRsolve(M,x,ep,alpha,rhs)
N = size(x,1); % Size of input data
r = size(rhs,2); % Number of right hand sides to consider
sol = zeros(M,r);

Dx = (1+(2*ep/alpha)^2)^.25*alpha*x;

phi = rbfphi(1:M,x,ep,alpha);

tmp = phi'*phi(:,[1,M-1:M]);
d = tmp(:,2)*sqrt(M-1)/sqrt(2) - ApplyTM(tmp(:,3));

% Start with the first columnn of the inverse
u0 = eye(M,1);
v0 = tmp(:,1);
s0 = sqrt(tmp(1,1));
q0 = phi(:,1)/s0;

sol(1,:) = (q0'*rhs)/s0;

% Apply the second column of the inverse
s2ok = sqrt(2);
c = -v0(2)/s2ok/s0^2;

u1 = [c*s2ok;1;zeros(M-2,1)];
v1 = s2ok*(ApplyTM(v0) + c*v0);
v1(M) = v1(M) - s2ok*(d'*u0);
s1 = sqrt(u1(1)*v1(1)+v1(2));
q1 = (s2ok*s0/s1)*(Dx.*q0 + c*q0);

sol(1:2,:) = sol(1:2,:) + u1(1:2)*((q1'*rhs)/s1);

% Iterate through the remaining columns
for k = 2:M-1
    tmp = ApplyTM(u1);
    c = -(v1'*tmp)/s1^2;
    s2ok = sqrt(2/k);
    bosk = (s1/s0)^2*sqrt((k-1)/k);
    
    u2 = s2ok*(tmp+c*u1) - bosk*u0;
    v2 = s2ok*(ApplyTM(v1)+c*v1) - bosk*v0;
    v2(M) = v2(M) - s2ok*(d'*u1);
    s2 = sqrt(u2'*v2);
    q2 = (s2ok*s1/s2)*(Dx.*q1+c*q1) - (bosk*s0/s2)*q0;
    
    sol = sol + u2*((q2'*rhs)/s2);

    u0 = u1; u1 = u2;
    v0 = v1; v1 = v2;
    q0 = q1; q1 = q2;
    s0 = s1; s1 = s2;
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
%
% The initial if block allows you to avoid extra computation if x is the
% same size as the last time this function was called
function TMx = ApplyTM(x)
persistent M Mvec;
if isempty(Mvec)
    M = size(x,1);
    Mvec = sqrt(1:M-1)'/sqrt(2);
elseif M~=size(x,1);
    M = size(x,1);
    Mvec = sqrt(1:M-1)'/sqrt(2);
end
TMx = [Mvec.*x(2:end);0] + [0;Mvec.*x(1:end-1)];
end
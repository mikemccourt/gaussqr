function [invU,Dvec,Q] = computeUinvHermite(x,M)
% function [invU,Dvec,Q] = computeUinvHermite(x,M)
% This is a test function to determine the viability of the fast QR
% technique on the Hermite polynomials.  If it works here, we'll consider
% it for the eigenfunction stuff.
% 
% Here we are computing the factorization P'*P = U'*D*U
% We return inv(U) and D, although D is stored as a vector
% This function assumes it has access to HermiteProd, which computes the
% Hermite polynomial of some necessary order
% The three term recurrence in place is 
%    p0(x) = alpha0,  p1(x) = (alpha1*x - beta1)*p0(x)
%    pk(x) = (alphak*x - betak)*p{k-1}(x)-gamma{k-1}p{k-2}(x)
% where alpha0 = 1, alphak = 2, betak = 0, gammak = 2(k-1)
%
% Inputs:  x - the points at which we are doing the regression
%              should be a column vector, only 1D for right now
%          M - the length of the regression
% Outputs: invU - the inverse of the U matrix, where P'P = U'DU is the LDL
%                 factorization
%          Dvec - a column vector of the D matrix diagonal
%          Q    - the matrix with orthonormal columns so that
%                 P = Q sqrt(D) U

N = size(x,1);
invU = zeros(M);
Dvec = zeros(M,1);
Q = zeros(N,M);

P = zeros(N,M); % Can we avoid computing this?
for k=1:M
    P(:,k) = HermiteProd(k-1,x);
end

% NOTE: Typo from LG paper, so d is negative of what was written
tmp = P'*P(:,[1,M-1:M]);
d = (M-1)*tmp(:,2) - ApplyTM(tmp(:,3));

% Start with the first columnns
u0 = eye(M,1);
v0 = tmp(:,1);

Dvec(1) = u0'*v0;
invU(:,1) = u0;

q0 = 1/sqrt(Dvec(1))*ones(N,1);
Q(:,1) = q0;

% Get the second columns
tmp = ApplyTM(u0,1);
c = -tmp'*v0/Dvec(1);
u1 = 2*(tmp+c*u0);
v1 = 2*(ApplyTM(v0) + c*v0);
v1(M) = v1(M) - 2*(d'*u0);

Dvec(2) = u1'*v1;
invU(:,2) = u1;

q1 = 2*sqrt(Dvec(1)/Dvec(2))*(x.*q0+c*q0);
Q(:,2) = q1;

% Go through the full iteration
% Is there less error by computing v separately at each turn,
% because then there is no inner product adding up small errors?
for k=3:M
    tmp = ApplyTM(u1,1);
    c = -tmp'*v1/Dvec(k-1);
    b = .25*Dvec(k-1)/Dvec(k-2);
    u2 = 2*(tmp+c*u1 - 2*b*u0);
    v2 = 2*(ApplyTM(v1)+c*v1 - 2*b*v0);
    v2(M) = v2(M) - 2*(d'*u1);
    
    Dvec(k) = u2'*v2;
    invU(:,k) = u2;
    
    q2 = 2*sqrt(Dvec(k-1)/Dvec(k))*(x.*q1+c*q1) - 4*b*sqrt(Dvec(k-2)/Dvec(k))*q0;
    Q(:,k) = q2;
    
    u0 = u1; u1 = u2;
    v0 = v1; v1 = v2;
    q0 = q1; q1 = q2;
end

end

% This is a private function that shouldn't be accessed outside this file
% It applies the matrix
%      TM = diag(1:M-1,-1)+diag(.5*ones(1,M-1),1)
% to a vector without forming that matrix
% M is the length of the vector
%
% If transp is passed at all, it applies the transpose of TM
function TMx = ApplyTM(x,transp)
M = size(x,1);
if exist('transp')
    TMx = [(1:M-1)'.*x(2:M);0] + [0;.5*x(1:M-1)];
else
    TMx = [.5*x(2:M);0] + [0;(1:M-1)'.*x(1:M-1)];
end
end
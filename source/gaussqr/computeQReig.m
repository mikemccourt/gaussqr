function [invU,Svec,Q] = computeQReig(M,x,ep,alpha,rhs)
% function [invU,Svec,Q] = computeQReig(M,x,ep,alpha,rhs)
%
% This computes the decomposition
%    P = Q * S * U
% where: P - gqr_phi(1:M,x,ep,alpha)
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
%                  gqr_phi(1:M,x,ep,alpha)*sol = rhs
%          without returning the factorization
%          Multiple rhs can be passed in matrix form
%
% Inputs: M - length of the eigenfunction summation
%             note that M<size(x,1)
%         x - the collocation points
%         ep - Gaussian shape parameter
%         alpha - GaussQR global scale parameter
%         rhs - right hand sides of Phi*sol = rhs
%
% Outputs: sol - column vectors such that Phi*sol = rhs
%
%%%%%
%
% Note that this may have less stability than the traditional QR
% factorization via Householder or whatever.  I still have some more
% stability analysis to do, but you should be aware that for larger M it
% could be an issue.

if nargin<4
    error('Insufficient arguments passed')
else
    [N,c] = size(x);
%     if M~=floor(M)
%         M = floor(M);
%         warning('Noninteger M passed, reset to M=%d',M)
%     end
%     
%     if M>=N
%         error('May only consider M<N for QR factorization, M=%g',M)
%     elseif M<1
%         error('Invalid value for M=%g',M)
%     elseif c>1
%         error('May only consider 1D problems, size(x,2)=%d',c)
%     end
    
    GQR = gqr_solveprep(1,x,ep,alpha,M);
    
% If the size of the problem is too small, the recurrence doesn't make
% sense, so direct evaluation is computed
    if nargin==4
        if M<4
            [Q,R] = qr(gqr_phi(GQR,x),0);
            Svec = diag(R);
            invU = triu(inv(diag(1./Svec)*R));
        else
            [invU,Svec,Q] = computeQRfactorization(GQR,x);
        end
    elseif nargin==5
        r = size(rhs,1);
        if N~=r
            error('dimension mismatch: size(x,1)=%d, size(rhs,1)=%d',N,r)
        elseif nargout>1
            error('Only one output returned when rhs is passed')
        end
        if M<4
            invU = linsolve(gqr_phi(GQR,x),rhs);
        else
            invU = computeQRsolve(GQR,x,rhs);
        end
    end
end
end


%%%%%%%%%%%
% This computes the actual factorization, and shouldn't be called directly
% from outside this file
function [invU,Svec,Q] = computeQRfactorization(GQR,x)
N = size(x,1);
M = GQR.Marr(end);
invU = zeros(M);
Svec = zeros(M,1);
Q = zeros(N,M);

Dx = (1+(2*GQR.ep/GQR.alpha)^2)^.25*GQR.alpha*x;

phi = gqr_phi(GQR,x);

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


%%%%%%%%%%%
% This private function computes the actual solution without storing the
% full factorization.  That allows for less storage.  Do not call this
% function from outside this file
function sol = computeQRsolve(GQR,x,rhs)
N = size(x,1); % Size of input data
M = GQR.Marr(end);
r = size(rhs,2); % Number of right hand sides to consider
sol = zeros(M,r);

Dx = (1+(2*GQR.ep/GQR.alpha)^2)^.25*GQR.alpha*x;

phi = gqr_phi(GQR,x);

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
    bosk = (s1/s0)^2*sqrt((k-1)/k); % this replaces b
    
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


%%%%%%%%%%%
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
% This shouldn't be a problem, because it will only ever be called on the
% column vectors of u and v, but I'm leaving a note in case I go senile.
%
% The if block allows you to avoid extra computation if x is the
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

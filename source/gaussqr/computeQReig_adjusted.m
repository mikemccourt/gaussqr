function [invU,Svec,Q] = computeQReig_adjusted(M,x,ep,alpha,rhs)
% function [invU,Svec,Q] = computeQReig_adjusted(M,x,ep,alpha,rhs)
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
    if M~=floor(M)
        M = floor(M);
        warning('Noninteger M passed, reset to M=%d',M)
    end
    
    if M>=N
        error('May only consider M<N for QR factorization, M=%g',M)
    elseif M<1
        error('Invalid value for M=%g',M)
    elseif c>1
        error('May only consider 1D problems, size(x,2)=%d',c)
    end
    
    [ep,alpha] = gqr_solveprep(1,x,ep,alpha);
    
% If the size of the problem is too small, the recurrence doesn't make
% sense, so direct evaluation is computed
    if nargin==4
        if M<4
            [Q,R] = qr(gqr_phi(gqr_formMarr(M),x,ep,alpha),0);
            Svec = diag(R);
            invU = triu(inv(diag(1./Svec)*R));
        else
            [invU,Svec,Q] = computeQRfactorization(M,x,ep,alpha);
        end
    elseif nargin==5
        r = size(rhs,1);
        if N~=r
            error('dimension mismatch: size(x,1)=%d, size(rhs,1)=%d',N,r)
        elseif nargout>1
            error('Only one output returned when rhs is passed')
        end
        if M<4
            invU = linsolve(gqr_phi(gqr_formMarr(M),x,ep,alpha),rhs);
        else
            invU = computeQRsolve(M,x,ep,alpha,rhs);
        end
    end
end
end


%%%%%%%%%%%
% This computes the actual factorization, and shouldn't be called directly
% from outside this file
function [invU,Svec,Q] = computeQRfactorization(M,x,ep,alpha)
N = size(x,1);
invU = zeros(M);
Svec = zeros(M,1);
Q = zeros(N,M);

Dx = (1+(2*ep/alpha)^2)^.25*alpha*x;

% Start with the first columnns
sq0 = gqr_phi(1,x,ep,alpha);
Svec(1) = norm(sq0);
Q(:,1) = sq0/Svec(1);

u0 = eye(M,1);
invU(:,1) = u0;

% Get the second columns
c = -(sq0'*(Dx.*sq0))/Svec(1)^2;
s2ok = sqrt(2);

sq1 = s2ok*(Dx.*sq0+c*sq0);
Svec(2) = norm(sq1);
Q(:,2) = sq1/Svec(2);

u1 = s2ok*(ApplyTM(u0)+c*u0);
invU(:,2) = u1;

% Iterate through the remaining columns
for k = 2:M-1
    c = -(sq1'*(Dx.*sq1))/Svec(k)^2;
    b = .5*(k-1)*(Svec(k)/Svec(k-1))^2;
    s2ok = sqrt(2/k);
    bosk = 2*b/sqrt(k^2-k);
    
    sq2 = s2ok*(Dx.*sq1+c*sq1) - bosk*sq0;
    Svec(k+1) = norm(sq2);
    Q(:,k+1) = sq2/Svec(k+1);
    
    u2 = s2ok*(ApplyTM(u1)+c*u1) - bosk*u0;
    invU(:,k+1) = u2;

    u0 = u1; u1 = u2;
    sq0 = sq1; sq1 = sq2;
end
end


%%%%%%%%%%%
% This private function computes the actual solution without storing the
% full factorization.  That allows for less storage.  Do not call this
% function from outside this file
function sol = computeQRsolve(M,x,ep,alpha,rhs)
N = size(x,1); % Size of input data
r = size(rhs,2); % Number of right hand sides to consider
sol = zeros(M,r);

Dx = (1+(2*ep/alpha)^2)^.25*alpha*x;

% Start with the first columnns
sq0 = gqr_phi(1,x,ep,alpha);
s0 = norm(sq0);
u0 = eye(M,1);

sol(1,:) = (sq0'*rhs)/s0^2;

% Apply the second column of the inverse
s2ok = sqrt(2);
c = -(sq0'*(Dx.*sq0))/s0^2;

sq1 = s2ok*(Dx.*sq0 + c*sq0);
s1 = norm(sq1);
u1 = [c*s2ok;1;zeros(M-2,1)];

sol(1:2,:) = sol(1:2,:) + u1(1:2)*((sq1'*rhs)/s1^2);

% Iterate through the remaining columns
for k = 2:M-1
    c = -(sq1'*(Dx.*sq1))/s1^2;
    s2ok = sqrt(2/k);
    bosk = (s1/s0)^2*sqrt((k-1)/k); % this replaces b
    
    u2 = s2ok*(ApplyTM(u1)+c*u1) - bosk*u0;
    sq2 = s2ok*(Dx.*sq1+c*sq1) - bosk*sq0;
    s2 = norm(sq2);
    
    sol = sol + u2*((sq2'*rhs)/s2^2);

    u0 = u1; u1 = u2;
    sq0 = sq1; sq1 = sq2;
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

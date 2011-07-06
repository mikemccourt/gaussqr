function X = ranksolve(U,VT,B)
% This solves a system (eye(n)+U*VT)*X=B
% where [n m]=size(U) and [m n]=size(VT)
% using the identity inv(eye(n)+U*VT)=eye(n)-U*inv(eye(m)+VT*U)*VT
% If m>.75n, the identity is not used.
[n m] = size(U);
[mV nV] = size(VT);
[nB mB] = size(B);

if n~=nV || m~=mV || n~=nB
    error('Unacceptable sizes in ranksolve, remember transpose')
end

if m>.75*n
    X = (eye(n)+U*VT)\B;
else
    X = B-U*((eye(m)+VT*U)\(VT*B));
end
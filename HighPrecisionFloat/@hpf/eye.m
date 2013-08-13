function mati = eye(varargin)
% hpf.eye - Identity matrix as HPF, of specified size and the default number of digits
% 
% hpf.eye(N) is the N-by-N identity matrix.
% 
%    hpf.eye(M,N) or hpf.eye([M,N]) is an M-by-N matrix with 1's on
%    the diagonal and zeros elsewhere.
% 
%    hpf.eye(SIZE(A)) is the same size as A.
% 
%    hpf.eye with no arguments is the scalar 1.

if (nargin == 0) || isempty(varargin{1})
  mati = hpf('1');
else
  % fill it with zeros, then overwrite the diagonal
  mati = repmat(hpf,varargin{:});
  
  % overwrite ones.
  mati(find(builtin('eye',varargin{:}))) = hpf('1');
end









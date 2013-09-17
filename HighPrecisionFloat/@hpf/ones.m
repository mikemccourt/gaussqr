function mat1 = ones(varargin)
% The matrix of all ones, of specified size and the default number of digits
%
% hpf.ones(N) is an N-by-N matrix of ones.
%  
% hpf.ones(M,N) or hpf.ones([M,N]) is an M-by-N matrix of ones.
%  
% hpf.ones(M,N,P,...) or hpf.ones([M N P ...]) is an M-by-N-by-P-by-... array of
% ones.
% 
% hpf.ones(SIZE(A)) is the same size as A and all ones.
% 
% hpf.ones with no arguments is the scalar hpf(1).

if (nargin == 0) || isempty(varargin{1})
  mat1 = hpf('1');
else
  mat1 = repmat(hpf('1'),varargin{:});
end


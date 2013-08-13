function mat0 = zeros(varargin)
% The matrix of all zeros, of specified size and the default number of digits
%
% hpf.zeros(N) is an N-by-N matrix of zeros.
%  
% hpf.zeros(M,N) or hpf.zeros([M,N]) is an M-by-N matrix of zeros.
%  
% hpf.zeros(M,N,P,...) or hpf.zeros([M N P ...]) is an M-by-N-by-P-by-... array of
% zeros.
% 
% hpf.zeros(SIZE(A)) is the same size as A and all zeros.
% 
% hpf.zeros with no arguments is the scalar 0.

if (nargin == 0) || isempty(varargin{1})
  mat0 = hpf;
else
  mat0 = repmat(hpf,varargin{:});
end



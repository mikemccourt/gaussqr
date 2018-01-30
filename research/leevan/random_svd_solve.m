function [x, logdet] = random_svd_solve(A, b, sval_thresh)

if nargin == 2
    sval_thresh = 0;
end
assert(sval_thresh >= 0)
assert(size(A, 1) == size(A, 2) && size(A, 1) == size(b, 1))

[U, S, V] = svd(A);
s_diag = diag(S);
good_svals = s_diag > sval_thresh;

x = V(:, good_svals) * bsxfun(@rdivide, U(:, good_svals)' * b, s_diag(good_svals));
if nargout == 2
    logdet = sum(log(s_diag(good_svals)));
end
end
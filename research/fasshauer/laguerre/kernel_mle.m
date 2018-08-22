% Could consider a Tikhonov parameter separately here, if we wanted to?
function mle = kernel_mle(K, y)
try
    R = chol(K);
catch exception
    if strcmp(exception.identifier, 'MATLAB:posdef')
        mle = 1e30;
        return
    end
    rethrow(exception)
end

N = length(K);
opts.UT = false;
opts.LT = true;
opts.TRANSA = true;
temp = linsolve(R, y, opts);

opts.UT = true;
opts.LT = false;
opts.TRANSA = false;
K_inv_y = linsolve(R, temp, opts);
mle = N * log(max(y' * K_inv_y, 1e-16)) + sum(log(max(diag(R), 1e-16)));
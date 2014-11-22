function val = StatFit1_CV_eval(ep,x_valid,y_valid,x_train,y_train,alpha)
% function val = StatFit1_CV_eval(ep,x_valid,y_valid,x_train,y_train,alpha)
% This function accepts sets of training and validation data and conducts
% the cross-validation using the HS-SVD
% Inputs: ep      - shape parameter for CV testing
%         x_valid - validation data points
%         y_valid - validation data values
%         x_train - training data points
%         y_train - training data values
%         alpha   - (optional) GQR scale parameter <default=1>
% Output: val     - cross-validation residual

switch nargin
    case 5
        alpha = 1;
    case 6
    otherwise
        error('Unnacceptable parameters, nargin=%d',nargin)
end

valid_num = length(x_valid);

val = 0;
for k=1:valid_num
    GQR = gqr_solve(x_train{k},y_train{k},ep,alpha);
    residual = norm(gqr_eval(GQR,x_valid{k}) - y_valid{k},2);
    val = val + residual^2;
end

val = sqrt(val);
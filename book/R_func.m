function R = R_func(a,x,z,vec_opt)
% function R = R_func(a,x,z,vec_opt)
% This function evaluates the boundary term from the GSS kernel derived by
% modifying the thin plate splines
% You can call this function to return a vector of R values or to return a
% matrix of R values, depending on if you just want to evaluate R or if you
% want to perform interpolation with it.
%
% Inputs: a - array of coefficients
%         x - matrix of 2D evaluation points
%         z - matrix of 2D kernel centers
%         vec_opt - evaluate R as a vector, not a matrix <default=0>
% Output: R - the boundary term evaluation
%
% R = R_func(a,x,z)
%         Assumes size(x,1) = N, size(z,1) = M
%         Returns a NxM matrix with x_i and z_j for the (i,j) entry
% R = R_func(a,x,z,1)
%         Assumes x and z both have N rows and 2 columns
%         Returns a Nx1 vector with x_i and z_i for the ith entry

if nargin<4
    vec_opt = 0;
end

if vec_opt==1
    % This assumes that you pass in 2D vectors of x and z points (which means
    % nx2 matrices) and an array of 3 values for a
    % It returns an nx1 vector R evaluated at each x and z combo
    R = 1/4*a(1)+12/29*(a(2)+a(3)) ...
        - 6/29*(a(2)*(x(:,1)+z(:,1)) + a(3)*(x(:,2)+z(:,2))) ...
        + 3/29*(a(2)*x(:,1).*z(:,1) + a(3)*x(:,2).*z(:,2));
else
    % This below is an attempt to create a matrix R(x,z) in the same way that
    % we need the matrix for interpolating
    r = size(x,1);
    c = size(z,1);
    x1 = repmat(x(:,1),1,c);
    x2 = repmat(x(:,2),1,c);
    z1 = repmat(z(:,1)',r,1);
    z2 = repmat(z(:,2)',r,1);
    R = 1/4*a(1)+12/29*(a(2)+a(3)) ...
        - 6/29*(a(2)*(x1+z1) + a(3)*(x2+z2)) ...
        + 3/29*(a(2)*x1.*z1 + a(3)*x2.*z2);
end
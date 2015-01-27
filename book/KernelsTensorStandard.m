% KernelsTensorStandard
% This is a straightforward example showing how a kernel can be evaluated
% in Matlab as the product of multiple kernels

% Define a series of 1D kernels which will be composed into a tensor
% kernel.  Note these are not defined in their radial form, even if some
% happen to be radial:
% The C2 Matern
KM2 = @(e,x,z) exp(-e*abs(bsxfun(@minus,x,z'))).*...
                 (1+e*abs(bsxfun(@minus,x,z')));
% The Inverse Multiquadric
KIM = @(e,x,z) 1./sqrt(1+(e*bsxfun(@minus,x,z')).^2);
% The Analytic Chebyshev with a = .5
% Note b is subbed for e, though they serve the same purpose
% Also note that 0<b<1
KCA = @(b,x,z) .5 + (1-b)* ...
         (b*(1-b^2) - 2*b*bsxfun(@plus,x.^2,z.^2') + (1+3*b^2)*x*z')./ ...
         ((1-b^2)^2 + 4*b*(b*bsxfun(@plus,x.^2,z.^2')-(1+b^2)*x*z'));
% The C0 Wendland kernel
KW0 = @(e,x,z) max(1-e*abs(bsxfun(@minus,x,z')),0);

% Organize these kernels into a cell array
Kcell = {KM2,KIM,KCA,KW0};

% Define the evaluation function for the kernel
% Because K is a 4D kernel, this function should accept
%          epvec - 1-by-4 vector
%          x - Neval-by-4 matrix
%          z - N-by-4 matrix
% The columns of each of the inputs get passed to the lower kernels for
% evaluation, and the results are combined in a product
Kf = @(e,x,z) prod(cell2mat(reshape( ...
         cellfun(@(K,e1,x1,z1) K(e1,x1,z1), ...
         Kcell,num2cell(e),num2cell(x,1), ...
         num2cell(z,1),'UniformOutput',0), ...
                            [1,1,length(e)])),3);
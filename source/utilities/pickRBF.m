function [rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF(type)
% At the moment, this function is only appropriate for 3D problems.  We'l
% make it more general eventually.
%
% We also need to fix this to throw an error if the user tries to call it
% with too many derivatives
% 
% Input:
% type   = function type (string: gaussian, imq, mq, linearmatern,
%                                 wendland_c2, wendland_c4, wendland_c6,
%                                 fundamental_3d)
% Output:
% rbf    = radial basis function;
% dxrbf  = camponent of the gradient of the RBF along x;
% dyrbf  = camponent of the gradient of the RBF along y;
% dzrbf  = camponent of the gradient of the RBF along z;
% Lrbf   = Laplacian of the RBF.
%
switch lower(type)
    case 'gaussian'
        rbf = @(e,r) exp(-(e*r).^2);
        dxrbf = @(e,r,dx) -dx*2*e^2.*exp(-(e*r).^2);
        dyrbf = @(e,r,dy) -dy*2*e^2.*exp(-(e*r).^2);
        dzrbf = @(e,r,dz) -dz*2*e^2.*exp(-(e*r).^2);
        Lrbf = @(e,r) 2*e^2*exp(-(e*r).^2).*(2*(e*r).^2-3);
    case 'imq'
        rbf = @(e,r) 1./sqrt(1+(e*r).^2);
        dxrbf = @(e,r,dx) -dx*e^2./(1+(e*r).^2).^(3/2);
        dyrbf = @(e,r,dy) -dy*e^2./(1+(e*r).^2).^(3/2);
        dzrbf = @(e,r,dz) -dz*e^2./(1+(e*r).^2).^(3/2);
        Lrbf = @(e,r) -3*e^2./(1+(e*r).^2).^(5/2);
    case 'mq'
        rbf = @(e,r) sqrt(1+(e*r).^2);
        dxrbf = @(e,r,dx) -dx*e^2./sqrt(1+(e*r).^2);
        dyrbf = @(e,r,dy) -dy*e^2./sqrt(1+(e*r).^2);
        dzrbf = @(e,r,dz) -dz*e^2./sqrt(1+(e*r).^2);
        Lrbf = @(e,r) e^2*(2*(e*r).^2+3)./(1+(e*r).^2).^(3/2);
    case 'linearmatern'
        rbf = @(e,r) exp(-(e*r)).*(1+e*r);
        dxrbf = @(e,r,dx) -dx*e^2.*exp(-(e*r));
        dyrbf = @(e,r,dy) -dy*e^2.*exp(-(e*r));
        dzrbf = @(e,r,dz) -dz*e^2.*exp(-(e*r));
        Lrbf = @(e,r) e^2*exp(-(e*r)).*(e*r-3);
    case 'wendland_c2'
        rbf = @(e,r) max(1-e*r,0).^4.*(4*e*r+1);
        dxrbf = @(e,r,dx) -dx*20*e^2.*max(1-e*r,0).^3;
        dyrbf = @(e,r,dy) -dy*20*e^2.*max(1-e*r,0).^3;
        dzrbf = @(e,r,dz) -dz*20*e^2.*max(1-e*r,0).^3;
        Lrbf = @(e,r) 60*e^2*(2*e*r-1).*max(1-e*r,0).^2;
    case 'wendland_c4'
        rbf = @(e,r) max(1-e*r,0).^6.*(35*(e*r).^2+18*(e*r)+3);
        dxrbf = @(e,r,dx) -dx*56*e^2.*(5*e*r+1).*max(1-e*r,0).^5;
        dyrbf = @(e,r,dy) -dy*56*e^2.*(5*e*r+1).*max(1-e*r,0).^5;
        dzrbf = @(e,r,dz) -dz*56*e^2.*(5*e*r+1).*max(1-e*r,0).^5;
        Lrbf = @(e,r) 56*e^2*(45*(e*r).^2-12*e*r-3).*max(1-e*r,0).^4;
    case 'wendland_c6'
        rbf = @(e,r) max(1-e*r,0).^8.*(32*(e*r).^3+25*(e*r).^2+8*e*r+1);
        dxrbf = @(e,r,dx) -dx*22*e^2.*(16*(e*r).^2+7*(e*r)+1).*max(1-e*r,0).^7;
        dyrbf = @(e,r,dy) -dy*22*e^2.*(16*(e*r).^2+7*(e*r)+1).*max(1-e*r,0).^7;
        dzrbf = @(e,r,dz) -dz*22*e^2.*(16*(e*r).^2+7*(e*r)+1).*max(1-e*r,0).^7;
        Lrbf = @(e,r) 22*e^2*(192*(e*r).^3-3*(e*r).^2-18*e*r-3).*max(1-e*r,0).^6;
    case 'fundamental_3d' % To use in the Method of Fundamental Solutions (MFS)
        rbf = @(r) 1./r;
        dxrbf = @(r,dx) -dx./r.^3;
        dyrbf = @(r,dy) -dy./r.^3;
        dzrbf = @(r,dz) -dz./r.^3;
        Lrbf = 0; % MFS doesn't require Laplacian!
    otherwise
        error('Function not recognized: check the input string')
end
end
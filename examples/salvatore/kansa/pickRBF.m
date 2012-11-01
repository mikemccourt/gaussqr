function [rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF(type)
% 
% Input:
% type   = function type (string: gaussian, imq, mq, linearmatern,
%                                 wendland_c2, wendland_c4, wendland_c6)
% Output:
% rbf    = radial basis function;
% dxrbf  = camponent of the gradient of the RBF along x;
% dyrbf  = camponent of the gradient of the RBF along y;
% dzrbf  = camponent of the gradient of the RBF along z;
% Lrbf   = Laplacian of the RBF
%

switch lower(type)
    case 'gaussian'
        rbf = @(e,r) exp(-(e.*r).^2);
        dxrbf = @(e,r,dx) -dx.*2*e.^2.*exp(-(e.*r).^2);
        dyrbf = @(e,r,dy) -dy.*2*e.^2.*exp(-(e.*r).^2);
        dzrbf = @(e,r,dz) -dz.*2*e.^2.*exp(-(e.*r).^2);
        Lrbf = @(e,r) -2*e.^2.*exp(-(e.*r).^2).*(2*(e.*r).^2+1);
    case 'imq' %
        rbf = @(e,r) 1./sqrt(1+(e.*r).^2);
        dxrbf = @(e,r,dx) -dx.*e.^2./(1+(e.*r).^2).^(3/2);
        dyrbf = @(e,r,dy) -dy.*e.^2./(1+(e.*r).^2).^(3/2);
        dzrbf = @(e,r,dz) -dz.*e.^2./(1+(e.*r).^2).^(3/2);
        Lrbf = @(e,r) -3*e.^2./(1+(e.*r).^2).^(5/2);
    case 'mq'
        rbf = @(e,r) sqrt(1+(e.*r).^2);
        dxrbf = @(e,r,dx) -dx.*e.^2./sqrt(1+(e.*r).^2);
        dyrbf = @(e,r,dy) -dy.*e.^2./sqrt(1+(e.*r).^2);
        dzrbf = @(e,r,dz) -dz.*e.^2./sqrt(1+(e.*r).^2);
        Lrbf = @(e,r) e.^2.*(2*(e.*r).^2+3)./(1+(e.*r).^2).^(3/2);
    case 'linearmaterne'
        rbf = @(e,r) exp(-(e.*r)).*(1+e.*r);
        dxrbf = @(e,r,dx) -dx.*e.^2.*exp(-(e.*r));
        dyrbf = @(e,r,dy) -dy.*e.^2.*exp(-(e.*r));
        dzrbf = @(e,r,dz) -dz.*e.^2.*exp(-(e.*r));
        Lrbf = @(e,r) e.^2.*exp(-(e.*r)).*(e.*r-3);
    case 'wendland_c2'
        rbf = @(e,r) (1-e.*r).^4.*(4*e.*r+1);
        dxrbf = @(e,r,dx) -dx.*20*e.^2.*max((1-e.*r).^3,0);
        dyrbf = @(e,r,dy) -dy.*20*e.^2.*max((1-e.*r).^3,0);
        dzrbf = @(e,r,dz) -dz.*20*e.^2.*max((1-e.*r).^3,0);
        Lrbf = @(e,r) 20.*e.^2.*(6*e.*r-3).*(1-e.*r).^2;
    case 'wendland_c4'
        rbf = @(e,r) (1-e*r).^6.*(35*(e.*r).^2+18*(e.*r)+3);
        dxrbf = @(e,r,dx) -dx.*56*e.^2.*(5*e.*r+1).*max((1-e.*r).^5,0);
        dyrbf = @(e,r,dy) -dy.*56*e.^2.*(5*e.*r+1).*max((1-e.*r).^5,0);
        dzrbf = @(e,r,dz) -dz.*56*e.^2.*(5*e.*r+1).*max((1-e.*r).^5,0);
        Lrbf = @(e,r) 56.*e.^2.*(45*(e.*r).^2-12*e.*r-3).*(1-e.*r).^4;
    case 'wendland_c6'
        rbf = @(e,r) (1-e.*r).^8.*(32*(e.*r).^3+25.*(e.*r).^2+8*e.*r+1);
        dxrbf = @(e,r,dx) -dx.*22*e.^2.*(16*(e.*r).^2+7*(e.*r)+1).*max((1-e.*r).^7,0);
        dyrbf = @(e,r,dy) -dy.*22*e.^2.*(16*(e.*r).^2+7*(e.*r)+1).*max((1-e.*r).^7,0);
        dzrbf = @(e,r,dz) -dz.*22*e.^2.*(16*(e.*r).^2+7*(e.*r)+1).*max((1-e.*r).^7,0);
        Lrbf = @(e,r) 22.*e.^2.*(192*(e.*r).^3-3*(e.*r).^2-18*e.*r-3).*(1-e.*r).^6;
    case 'fund'
        rbf = @(e,r) 1./r;
        dxrbf = @(e,r,dx) -dx./r.^3;
        dyrbf = @(e,r,dy) -dy./r.^3;
        dzrbf = @(e,r,dz) -dz./r.^3;
        Lrbf = @(e,r) 0*r;
end

end
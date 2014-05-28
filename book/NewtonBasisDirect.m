% NewtonBasisDirect
% Computes and plots Newton basis functions for 2D RBF interpolation
% Calls on: DistanceMatrix
  rbf = @(e,r) exp(-(e*r).^2); ep = 3;
  %rbf = @(e,r) sqrt(1+(e*r).^2); ep = 5;  % MQ RBF
  %rbf = @(epsilon,r) exp(-epsilon*r);  ep=1;  % basic Matern
  N = 25; gridtype = 'u';
  dsites = CreatePoints(N,2,gridtype);
  ctrs = dsites;
  neval = 40; M = neval^2;
  epoints = CreatePoints(M,2,'u');
  DM_data = DistanceMatrix(dsites,ctrs);
  IM = rbf(ep,DM_data);
  % Better to use regular evaluation matrix. Then
  DM_B = DistanceMatrix(epoints,ctrs);
  B = rbf(ep,DM_B);
  % Find all Newton functions directly without iteration (works only for
  % positive definite kernels)
  L = chol(IM,'lower');  % produces lower triangular matrix L
  beta = L*diag(diag(L)); % scale columns
  newtonbasis = (beta\B')';
  %test = (beta\IM')'; % should be lower triangular matrix
  figure
  xe = reshape(epoints(:,1),neval,neval);
  ye = reshape(epoints(:,2),neval,neval);
  NFplot = surf(xe,ye,reshape(newtonbasis(:,ceil(N/2)),neval,neval));
  %ctrs(ceil(N/2),:)
  set(NFplot,'FaceColor','interp','EdgeColor','none')
  colormap autumn; view([145 45]); camlight; lighting gouraud
  figure
  NFplot = surf(xe,ye,reshape(newtonbasis(:,1),neval,neval));
  %ctrs(1,:)
  set(NFplot,'FaceColor','interp','EdgeColor','none')
  colormap autumn; view([145 45]); camlight; lighting gouraud
  figure
  NFplot = surf(xe,ye,reshape(newtonbasis(:,ceil(sqrt(N)/2)),neval,neval));
  %ctrs(ceil(sqrt(N)/2),:)
  set(NFplot,'FaceColor','interp','EdgeColor','none')
  colormap autumn; view([145 45]); camlight; lighting gouraud
  
%   figure
%   idx=randi(N,1);
%   NFplot = surf(xe,ye,reshape(newtonbasis(:,idx),neval,neval));
%   %ctrs(idx,:)
%   set(NFplot,'FaceColor','interp','EdgeColor','none')
%   colormap autumn; view([145 45]); camlight; lighting gouraud
%   
  

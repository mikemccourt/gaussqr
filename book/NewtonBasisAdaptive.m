% NewtonBasisAdaptive
% Computes and plots Newton basis functions for 2D RBF interpolation
% Calls on: DistanceMatrix
  rbf = @(e,r) exp(-(e*r).^2); ep = 5;
  %rbf = @(epsilon,r) exp(-epsilon*r);  ep=1;  % basic Matern
  N = 25; gridtype = 'h';
  dsites = CreatePoints(N,2,gridtype);
  ctrs = dsites;
  neval = 40; M = neval^2;
  epoints = CreatePoints(M,2,'u');
  DM = DistanceMatrix(dsites,ctrs);
  Kinterp = rbf(ep,DM);
  % Better to use regular evaluation matrix. Then
  DM = DistanceMatrix(epoints,ctrs);
  Keval = rbf(ep,DM);
  newtonbasis = zeros(N,N);
  newtonbasisplot = zeros(M,N);
  % This will be the matrix containing the Newton basis, columnwise,
  % as functions on the evaluation grid.
  iter = 1; % iteration counter
  z = diag(Kinterp);
  [zmax, zind] = max(abs(z));
  newtonbasis(:,1) = Kinterp(:,zind)/sqrt(z(zind));
  newtonbasisplot(:,1) = Keval(:,zind)/sqrt(z(zind));
  w = newtonbasis(:,1).^2;
  % Now the main iteration. Each step uses a new point 
  % and a computes new basis function.
  for iter=2:N
    iter
    [zmax, zind] = max(z-w);
    y = newtonbasis(zind,1:iter-1)';
    u = Kinterp(:,zind) - newtonbasis(:,1:iter-1)*y;
    powerfuniter = sqrt(z(zind)-w(zind));
    newtonbasis(:,iter) = u/powerfuniter;
    w = w + newtonbasis(:,iter).^2;
    u = Keval(:,zind) - newtonbasisplot(:,1:iter-1)*y;
    newtonbasisplot(:,iter) = u/powerfuniter;
  end

%  figure
  xe = reshape(epoints(:,1),neval,neval);
  ye = reshape(epoints(:,2),neval,neval);
  CFplot = surf(xe,ye,reshape(newtonbasisplot(:,ceil(N/2)),neval,neval));
  set(CFplot,'FaceColor','interp','EdgeColor','none')
  colormap autumn; view([145 45]); camlight; lighting gouraud
  figure
  CFplot = surf(xe,ye,reshape(newtonbasisplot(:,1),neval,neval));
  set(CFplot,'FaceColor','interp','EdgeColor','none')
  colormap autumn; view([145 45]); camlight; lighting gouraud
  figure
  CFplot = surf(xe,ye,reshape(newtonbasisplot(:,ceil(sqrt(N)/2)),neval,neval));
  set(CFplot,'FaceColor','interp','EdgeColor','none')
  colormap autumn; view([145 45]); camlight; lighting gouraud
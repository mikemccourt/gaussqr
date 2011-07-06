function srat = rbfsvdF(Marr,x,ep,a,b,c)
rbfsvd = svd(rbfphi(Marr,x,ep,a,b,c));
srat = log10(rbfsvd(1)/rbfsvd(end)+eps);
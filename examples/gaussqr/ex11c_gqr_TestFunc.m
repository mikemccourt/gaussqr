function errval = ex11c_gqr_TestFunc(ep,epstruct)
    Marr = epstruct.Marr;
    alpha = epstruct.alpha;
    AMFpts = epstruct.AMFpts;
    BMFpts = epstruct.BMFpts;
    FmatMF = epstruct.FmatMF;
    FAcou = epstruct.FAcou;
    FAifa = epstruct.FAifa;
    FBcou = epstruct.FBcou;
    FBifa = epstruct.FBifa;
    Frhs = epstruct.Frhs;
    restart = epstruct.restart;
    Fsol = epstruct.Fsol;
    Fx = epstruct.Fx;
    try
        FmatFD = epstruct.FmatFD;
    catch
        FmatFD = [];
    end

    phiA = rbfphi(Marr,AMFpts,ep,alpha);
    phiAx = rbfphi(Marr,Fx(FBifa,:),ep,alpha,[1 0]);
    phiB = rbfphi(Marr,BMFpts,ep,alpha);
    phiBx = rbfphi(Marr,Fx(FBifa,:),ep,alpha,[1 0]);

    FmatMF(FBifa,:) = zeros(size(FmatMF(FBifa,:)));
    FmatMF(FBifa,[FAcou,FAifa]) = phiAx/phiA;
    FmatMF(FBifa,[FBcou,FBifa]) = -phiBx/phiB;
    Frhs(FBifa) = zeros(size(FBifa));

    [Funew,cnvg] = gmres(FmatMF,Frhs,restart,[],[],FmatFD);
    errval = errcompute(Funew,Fsol);
end

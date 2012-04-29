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

    phiA = gqr_phi(Marr,AMFpts,ep,alpha);
    phiAx = gqr_phi(Marr,Fx(FBifa,:),ep,alpha,[1 0]);
    phiB = gqr_phi(Marr,BMFpts,ep,alpha);
    phiBx = gqr_phi(Marr,Fx(FBifa,:),ep,alpha,[1 0]);

    FmatMF(FBifa,:) = zeros(size(FmatMF(FBifa,:)));
    FmatMF(FBifa,[FAcou,FAifa]) = phiAx/phiA;
    FmatMF(FBifa,[FBcou,FBifa]) = -phiBx/phiB;
    Frhs(FBifa) = zeros(size(FBifa));

    [Funew,cnvg] = gmres(FmatMF,Frhs,restart,[],[],FmatFD);
    errval = errcompute(Funew,Fsol);
end

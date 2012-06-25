function errval = ex11c_gqr_TestFunc(ep,epstruct)
% Note that in this function, epstruct can get passed to gqr_phi
% This is because gqr_phi only looks for the fields Marr, alpha and ep,
%  which epstruct has after I add ep below
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
    epstruct.ep = ep;

    phiA = gqr_phi(epstruct,AMFpts);
    phiAx = gqr_phi(epstruct,Fx(FBifa,:),[1 0]);
    phiB = gqr_phi(epstruct,BMFpts);
    phiBx = gqr_phi(epstruct,Fx(FBifa,:),[1 0]);

    FmatMF(FBifa,:) = zeros(size(FmatMF(FBifa,:)));
    FmatMF(FBifa,[FAcou,FAifa]) = phiAx/phiA;
    FmatMF(FBifa,[FBcou,FBifa]) = -phiBx/phiB;
    Frhs(FBifa) = zeros(size(FBifa));

    [Funew,cnvg] = gmres(FmatMF,Frhs,restart,[],[],FmatFD);
    errval = errcompute(Funew,Fsol);
end

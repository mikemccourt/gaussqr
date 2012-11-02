function rel_error = formula_error(R, pnts, srcpnt, dipmom, sigma, maxdegree)

analytic = HomSpherePotential(R(2), sigma, srcpnt, dipmom, pnts);
semianalytic = MultiSpherePotential(R,[sigma, sigma],srcpnt,dipmom,pnts,maxdegree);
rel_error = norm(analytic - semianalytic)/norm(analytic);

end

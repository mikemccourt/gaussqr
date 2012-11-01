% clear all

srcpnt = zeros(50,3);
srcpnt(:,3) = linspace(0.0001,0.0699,50);
points = SphereRegGoldPoints(1500, 0.1);
sigma = 0.2;
R = [0.07 0.1];
maxdegree = [30 50 70 90 100 110];

for j = 1: size(maxdegree,2);
    % Error versus dipole radial position along z, radially oriented dipole
    dipmom = 2.7e-12 / sqrt(3) * [1, 1, 1];
    rel_err1 = zeros(size(srcpnt,1),1);
    for i = 1:size(srcpnt,1)
        rel_err1(i) = formula_error_analytic_dip(R, points, srcpnt(i,:), dipmom, sigma, maxdegree(j));
    end
    
    % Error versus dipole radial position along z, x oriented dipole
    dipmom = 2.7e-12 * [1, 0, 0];
    rel_err2 = zeros(size(srcpnt,1),1);
    for i = 1:size(srcpnt,1)
        rel_err2(i) = formula_error_analytic_dip(R, points, srcpnt(i,:), dipmom, sigma, maxdegree(j));
    end
    
    % Error versus dipole position, dipole in a plane at 45 degrees with respect
    % to the xy plane
    dipmom = 2.7e-12 * [sqrt(2)/2, sqrt(2)/2, 0];
    rel_err3 = zeros(size(srcpnt,1),1);
    for i = 1:size(srcpnt,1)
        rel_err3(i) = formula_error_analytic_dip(R, points, srcpnt(i,:), dipmom, sigma, maxdegree(j));
    end
    
    createfigure_analytic_dip(srcpnt(:,3),[rel_err1, rel_err2, rel_err3]);
    j
end
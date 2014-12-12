function ctrs = centers( pnts, coeff )
% CENTERS generates the distribution of centers "ctrs" by 
% inflating/deflating the collocation points cloud "pnts" by a coefficient
% "coeff"

N = size(pnts,1);

% Calculate the centroid of the point cloud
centroid = sum(pnts,1); 
centroid = centroid/N;

% Reference the points to the centroid
ctrs(:,1) = pnts(:,1) - centroid(1);
ctrs(:,2) = pnts(:,2) - centroid(2);
ctrs(:,3) = pnts(:,3) - centroid(3);

% Find centers by inflation/deflation of the point cloud
ctrs = ctrs*coeff;

% Restore the original reference
ctrs(:,1) = ctrs(:,1) + centroid(1);
ctrs(:,2) = ctrs(:,2) + centroid(2);
ctrs(:,3) = ctrs(:,3) + centroid(3);
end


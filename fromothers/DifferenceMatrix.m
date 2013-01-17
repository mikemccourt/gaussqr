% DM = DifferenceMatrix(datacoord,centercoord)
% Borrowed from Fasshauer's "Meshfree Approximation Methods in Matlab"
%
% Forms the difference matrix of two sets of points in R
% (some fixed coordinate of point in R^s), i.e.,
% DM(j,k) = datacoord_j - centercoord_k .
  function DM = DifferenceMatrix(datacoord,centercoord)
  % The ndgrid command produces two MxN matrices:
  %   dr, consisting of N identical columns
  %       (each containing the M data sites)
  %   cc, consisting of M identical rows
  %       (each containing the N centers)
  
  % Test to allow for user to pass [] and return []
  if (length(datacoord)==0 || length(centercoord)==0)
      DM = [];
  else
      [dr,cc] = ndgrid(datacoord(:),centercoord(:));
      DM = dr-cc;
  end
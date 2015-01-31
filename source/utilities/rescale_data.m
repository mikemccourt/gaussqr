function [x_clean,y_clean,shift,scale,bad_ind] = rescale_data(x,y,keepunique)
% function [x_clean,y_clean,shift,scale,bad_ind] = rescale_data(x,y,keepunique)
% This function accepts in x data and rescales it to [-1,1]^d
% It also removes any NaN or Inf data
% Input:   x - real valued matrix size N-by-d
%          y - associated data values at locations x
%          keepunique - run unique on the x locations first
%                       (optional, default=0)
% Outputs: x_clean - x data, cleaned and rescaled into [-1,1]^d
%          y_clean - y data with bad values removed
%          shift - shift required to recover the data
%          scale - scale required to recover the data
%          bad_ind - list of the removed indices in the original data
% 
% To recover the original data points (minus any NaN or Inf),
%    x = scale*(x_scaled+1)/2 + shift;
%
% You can call this with just x locations and it will scale just the
% locations to [-1,1]^d:
% function x_in_minus1_plus1 = rescale_data(x)
%
% DEVELOPER'S NOTE: I'm not sure if it makes more sense to move the unique
% test to after checking for NaN and inf - probably the same either way.

if not(exist('y','var'))
    if nargout>1
        error('No y value was passed, no y value can be returned')
    else
        y = ones(size(x,1),1);
    end
end

if not(exist('keepunique','var'))
    keepunique = 0;
end

if keepunique
    [xunique,uniqueind] = unique(x,'rows');
    x = xunique;
    y = y(uniqueind);
end

% Find any troubling values and remove them
good_ind = not(or(isnan(sum([x,y],2)),isinf(sum([x,y],2))));
x = x(good_ind,:);
y_clean = y(good_ind);

N = size(x,1);

shift = min(x);
scale = max(x) - min(x);
x_clean = 2*(x - ones(N,1)*shift)./(ones(N,1)*scale) - 1;

bad_ind = find(not(good_ind));
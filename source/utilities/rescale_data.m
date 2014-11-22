function [x_clean,y_clean,shift,scale,bad_ind] = rescale_data(x,y);
% function [x_clean,y_clean,shift,scale,bad_ind] = rescale_data(x,y)
% This function accepts in x data and rescales it to [-1,1]^d
% It also removes any NaN or Inf data
% Input:   x - real valued matrix size N-by-d
% Outputs: x_clean - x data, cleaned and rescaled into [-1,1]^d
%          y_clean - y data with bad values removed
%          shift - shift required to recover the data
%          scale - scale required to recover the data
%          bad_ind - list of the removed indices in the original data
% 
% To recover the original data points (minus any NaN or Inf),
%    x = scale*(x_scaled+1)/2 + shift;

% Find any troubling values and remove them
good_ind = not(or(isnan(sum([x,y],2)),isinf(sum([x,y],2))));
x = x(good_ind,:);
y_clean = y(good_ind);

N = size(x,1);

shift = min(x);
scale = max(x) - min(x);
x_clean = 2*(x - ones(N,1)*shift)./(ones(N,1)*scale) - 1;

bad_ind = find(not(good_ind));
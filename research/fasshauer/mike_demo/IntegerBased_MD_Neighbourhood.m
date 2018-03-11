%-------------------------------------------------------------------------%
%
% File: IntegerBased_MD_Neighbourhood(dsites,idx_ds,index,q,M,t)
%
% Goal: script that founds the neighbouring blocks
%
% Inputs:  dsites:       NXM matrix representing a set of N data sites
%          idx_ds:       the integer-based data structure
%          index:        indices of points lying in the k-th block
%          q:            number of blocks in one direction
%          M:            space dimension
%          t:            it is related to the number of neighbouring 
%                        blocks, the procedure searches in the k-th block
%                        and in its (3 + 2^t)^M - 1 neighboring blocks
%
% Outputs:  dxx:         points lying in the k-th neighbourhood
%           dx:          indices of points lying in the k-th neighbourhood
%
%-------------------------------------------------------------------------%
function [dxx dx] = IntegerBased_MD_Neighbourhood(dsites,idx_ds,index,...
    q,M,t)
neigh = []; k = M-1; k1 = 1:t; % Initialize
% Find neighbouring blocks
while k  > 0
     neigh = [index+k1.*q.^k,index-k1.*q.^k];
    if k - 1 > 0
        neigh = [neigh,neigh+q.^(k-1),neigh-q.^(k-1)];
    end
    k = k - 1;
end
k2 = 1:t; k3 = k2;
for i = 1:length(neigh)
    neighplus(k2) = neigh(i) + k3;
    neighminus(k2) = neigh(i) - k3;
    k2 = k2 + t;
end
neigh = [neigh,index+k1,index-k1,neighplus,neighminus];
% Reduce the number of neighbouring blocks for border blocks
j = find(neigh > 0 & neigh <= q^M);
index = [index; neigh(j)'];
dxx = []; dx = []; 
for p = 1:length(index)
    dxx = [dxx;dsites(idx_ds{index(p)},:)];
    dx = [dx;idx_ds{index(p)}];
end
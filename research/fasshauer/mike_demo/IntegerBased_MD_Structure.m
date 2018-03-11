%-------------------------------------------------------------------------%
%
% File: IntegerBased_MD_Structure(dsites,q,puradius,M)
%
% Goal: find the data sites located in each of the q^M blocks
%
% Inputs: dsites:     NXM matrix representing a set of N data sites
%         q:          number of blocks in one direction
%         puradius:   radius of PU subdomains
%         M:          space dimension
%
% Outputs: idx_dsites_k: multiarray containing the indices of the data 
%                        points located in k-th block
%
%-------------------------------------------------------------------------%
function [idx_dsites_k] = IntegerBased_MD_Structure(dsites,q,puradius,M)
N = size(dsites,1); idx_dsites_k = cell(q^M,1); k = 1:M-1;
for i = 1:N
    idx = ceil(dsites(i,:)./puradius);
    idx(idx == 0) = 1;
    index = sum((idx(k)-1).*q.^(M-k)) + idx(end);
    idx_dsites_k{index} = [idx_dsites_k{index}; i];
end
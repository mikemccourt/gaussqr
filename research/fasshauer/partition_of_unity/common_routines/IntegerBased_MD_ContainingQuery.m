%-------------------------------------------------------------------------%
%
% File: IntegerBased_MD_ContainingQuery(puctr,q,puradius,M)
%
% Goal: script that given a subdomain centre returns the index of
%       the square block containing the subdomain centre
%
% Inputs: puctr:       subdomain centre
%         q:           number of blocks in one direction
%         puradius:    radius of the PU subdomains
%         M:           space dimension
%
% Outputs: index: the index of the block containing the subdomain
%                 centre
%
%-------------------------------------------------------------------------%
function [index] = IntegerBased_MD_ContainingQuery(puctr,q,puradius,M)
idx = ceil(puctr./puradius); k = 1:M-1;
idx(idx == 0) = 1;
index = sum((idx(k)-1).*q.^(M-k)) + idx(end);
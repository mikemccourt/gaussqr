%-------------------------------------------------------------------------%
%
% File: IntegerBased_MD_RangeSearchAniso(puctr,puradius,dsites,index)
%
% Goal: find the data sites located in a given subdomain 
%
% Inputs: puctr:     subdomain centre
%         puradius:  radius of PU subdomains
%         dsites:    NXM matrix representing a set of N data sites
%         index:     vector containing the indices of the data points 
%                    located in the k-th block (the block containing the 
%                    subdomain centre) and in the neighbouring blocks
%  
% Outputs: idx:  vector containing the indices of the data points located
%                in a given PU subdomain
%
%-------------------------------------------------------------------------%
function [idx] = IntegerBased_MD_RangeSearchAniso(puctr,puradius,...
    dsites,index)
[N,M] = size(dsites); ts = []; idx = []; % Initialize
% Compute distances between the data sites and the centre
k = 1;
for i = 1:N
    tot = 0;
    for d = 1:M
        tot = tot + ((dsites(i,d)-puctr(d))./puradius(d)).^2;
    end
    ts(i) = tot <= 1;
    %ts(i) = norm((dsites(i,:)-puctr)./puradius).^2 <= 1;
    if ts(i) > 0
        idx(k) = index(i);
        k = k + 1;
    end
end
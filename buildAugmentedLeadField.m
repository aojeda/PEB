function [H, Delta, blocks, indG, indV] = buildAugmentedLeadField(hm, M, chanlocs)
% Augments the lead field matrix with artifact scalp maps projections.
% 
% Inputs:
%   hm: is a head model object
%   M: is a dictionary of artifact scalp map projections (each column is a scalp map)
%   chanlocs: is the channel locations structure that corresponds to each column of M 
%
% Outputs:
%   H: [L A] is the lead field matrix augmented with a dictionary of artifact scalp maps
%   indG are the indices of the brain sources
%   indV are the indices of the artifact sources

xyz =  [[chanlocs.X]'  [chanlocs.Y]'  [chanlocs.Z]'];
labels  = {chanlocs.labels};
[~, ~, loc2] = intersect(hm.labels, labels,'stable');
[Ny,Ng] = size(hm.K);
Nic = size(M,2);
Nroi = length(hm.atlas.label);
A = zeros(Ny,Nic);
for ic=1:Nic
    F = scatteredInterpolant(xyz(loc2,:),M(loc2,ic));
    A(:,ic) = F(hm.channelSpace);
end
norm_K = norm(hm.K);
L = hm.K/norm_K;
Delta = hm.L/norm_K;
H = [L A];
H = bsxfun(@rdivide,H,sqrt(sum(H.^2)));
blocks = hm.indices4Structure(hm.atlas.label);
blocks = [[blocks;false(Nic,Nroi)] [false(Ng,Nic);diag(true(Nic,1))]];
indG = (1:Ng)';
indV = Ng+(1:Nic)';
end
function [clonesizes_sampled] = subsampling_clones_overtime(nsamples,sizesamples,clonesizes,filter1,timepoints)
%% Repeated sampling (random permutation) of simulated clones in subsets:
% It generates subsets with a limited number of randomly-chosen clones from
% the pool of simulated clones by random permutation (with replacement).
% Parameters nsamples and sizesamples specify the total number of subsets
% and the number of clones contained in each, respectively.

% from Colom et al, 2020

%% Input:
% nsamples: number of subsets the data is distributed in.
% sizesamples: number of random clones assigned per subset.
% clonesizes: cell array {1,timepoints}(:,1) with input clone size data at given time points
% filter1: threshold for clone sampling (only clones with size >= filter1 are sampled)
% timepoints: time points when clone sizes were calculated (expressed in weeks)

%% Output:
% clonesizes_sampled: cell array {subsets,1}{1,timepoints}(:,1) of clone sizes in each subset at each specified time point

%% Sampling (random permutation) of experimental clones:
clonesizes_sampled = {};

for luptime = 1:size(timepoints,2)
    % Restrict sampling to clones with a number of cells >= filter1
    loc_prolif = find(clonesizes{1,luptime}(:,1)>=filter1);
    rnd_pickpos = [];
    % Subset making:
    for rnd_loop = 1:nsamples
        for rnd_loop2 = 1:sizesamples(luptime)
            % random clone sampling:
            rnd_pickpos = loc_prolif(randperm(size(loc_prolif,1),1),1);
            clonesizes_sampled{rnd_loop,1}{1,luptime}(rnd_loop2,1) = clonesizes{1,luptime}(rnd_pickpos,1);
        end
    end
end

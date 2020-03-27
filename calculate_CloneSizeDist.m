function [CUM_nfreq,CUM_nfreq_centre,CUM_nfreq_ci95up,CUM_nfreq_ci95dn,avgCloneSize] = calculate_CloneSizeDist(timepoints,time2calc,cloneSizes,showCI,sampling)
%% Calculates the clone size distribution over the time series.
% The time course of the normalized clone size frequencies is obtained from
% the clone sizes at given time points.

% from Colom et al, 2020

%% Input:
% timepoints: vector of time points in the original data (expressed in weeks)
% time2calc: vector of time points of interest to retrieve data from (expressed in weeks)
% cloneSizes: matrix of size [m,n,q] containing clone sizes (m = No. of replicates, n = No. of time point, q = all clones)
% showCI: calculate plausible intervals on model outcome (based on experimental sampling error) ( 0=NO | 1=YES )
% sampling: structure containing sampling properties for CI estimation
    % struct{NSubsets==20, NClones==[6000 6000 6000]}
        % NSubsets: number of subsets simulation data are distributed to (random sampling with replacement)
        % NClones: total number of clones contained per time point in each subset

%% Output:
% CUM_nfreq: cell array of normalized cumulative clone size frequencies at given time points
% CUM_nfreq_centre: mid of 95 bound on normalized cumulative clone size frequencies (sampling error)
% CUM_nfreq_ci95up: upper 95% bound on normalized cumulative clone size frequencies (sampling error)
% CUM_nfreq_ci95dn: bottom 95% bound on normalized cumulative clone size frequencies (sampling error)
% avgCloneSize: average clone size over time

%% Example:
% timepoints = [1:3];
% time2calc = [2,3];
% cloneSizes(1,1,:) = [1 1 1]'; cloneSizes(1,2,:) = [5 5 5]'; cloneSizes(1,3,:) = [10 10 10]';
% showCI = 0;
% sampling = struct('NSubsets',20,'NClones',[6000 6000 6000]);
% [CUM_nfreq,CUM_nfreq_ci95up,CUM_nfreq_ci95dn,avgCloneSize] = calculate_CloneSizeDist(timepoints,time2calc,cloneSizes,showCI,sampling);

%% Default parameter values:
% Fill in unset optional values:
if (nargin < 3)
    showCI = 0;
    sampling = struct('NSubsets',20,'NClones',[11552 15865 4152 2474 3485]);
end

CUM_nfreq = {};
CUM_nfreq_centre = {};
CUM_nfreq_ci95up = {};
CUM_nfreq_ci95dn = {};

%% Calculate clone size distributions:
for aja = 1:length(time2calc)
    myTime = find(timepoints >= time2calc(aja),1); % position in the time vector
    cloneSizes_2plot{1,aja}(:,1) = cloneSizes(1,myTime,:);
    avgCloneSize(aja) = mean(cloneSizes_2plot{1,aja}(find(cloneSizes_2plot{1,aja}~=0)));
    nfreq_all{1,aja} = histc(cloneSizes_2plot{1,aja},[0:1:max(cloneSizes_2plot{1,aja})]);
    nfreq_all_rel{1,aja} = nfreq_all{1,aja}./sum(nfreq_all{1,aja}(2:end,1)); % only persistent clones considered
    CUM_nfreq{1,aja} = cumsum(nfreq_all_rel{1,aja}(2:end,1),'reverse');
end

%% Calculate plausible intervals between simulation runs (by clonal subsampling):
if showCI == 1
    % Subsampling with replacement: Nclones are sampled at each time point (to ~fit experimental data)
    cloneSizes_2plot_sampled = subsampling_clones_overtime(sampling.NSubsets,sampling.NClones,cloneSizes_2plot,1,time2calc);
    % Calculate clone size distributions from subsamples:
    for aja = 1:length(time2calc)
        nfreq_samples_rel = [];
        for eje = 1:sampling.NSubsets
            nfreq_sample = histc(cloneSizes_2plot_sampled{eje,1}{1,aja},[0:1:max(cloneSizes_2plot{1,aja})]);
            nfreq_samples_rel(1:size(nfreq_sample,1),eje) = nfreq_sample./sum(nfreq_sample(2:end,1)); % only persistent clones considered
        end
        nfreq_samples_all_rel{1,aja} = nfreq_samples_rel;
        CUM_nfreq_samples_all_rel{1,aja} = cumsum(nfreq_samples_all_rel{1,aja}(2:end,:),'reverse');        
        % Calculate 95% plausible intervals on the clone size freqs.:
        CUM_nfreq_ci95up{1,aja} = quantile(CUM_nfreq_samples_all_rel{1,aja},0.975,2);
        CUM_nfreq_ci95dn{1,aja} = quantile(CUM_nfreq_samples_all_rel{1,aja},0.025,2);
        CUM_nfreq_centre{1,aja} = quantile(CUM_nfreq_samples_all_rel{1,aja},0.5,2);
        for pos = 1:size(CUM_nfreq_centre{1,aja},1)
            if CUM_nfreq_ci95up{1,aja}(pos,1)<=0; CUM_nfreq_ci95up{1,aja}(pos,1)=1E-8; end
            if CUM_nfreq_ci95dn{1,aja}(pos,1)<=0; CUM_nfreq_ci95dn{1,aja}(pos,1)=1E-8; end
        end
    end
end

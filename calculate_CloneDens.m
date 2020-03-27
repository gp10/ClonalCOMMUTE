function [NClones,NClones_ci95up,NClones_ci95dn] = calculate_CloneDens(timepoints,cloneSizes,showCI,sampling)
%% Calculates the number of surviving clones over the time series.
% The time course of the number of surviving clones is obtained from the
% clone sizes at given time points.

% from Colom et al, 2020

%% Input:
% timepoints: vector of time points (expressed in weeks)
% cloneSizes: matrix of size [m,n,q] containing clone sizes (m = No. of replicates, n = No. of time point, q = all clones)
% showCI: calculate plausible intervals on model outcome (based on experimental sampling error) ( 0=NO | 1=YES )
% sampling: structure containing sampling properties for CI estimation
    % struct{NSubsets==20, NClones==6000}
        % NSubsets: number of subsets simulation data are distributed to (random sampling with replacement)
        % NClones: total number of clones contained in each subset

%% Output:
% NClones: vector with the number of surviving clones at different time points
% NClones_ci95up: upper 95% bound on the number of surviving clones (sampling error)
% NClones_ci95dn: bottom 95% bound on the number of surviving clones (sampling error)

%% Example:
% timepoints = [1:3];
% cloneSizes(1,1,:) = [1 1 1]'; cloneSizes(1,2,:) = [5 5 5]'; cloneSizes(1,3,:) = [10 10 10]';
% showCI = 0;
% sampling = struct('NSubsets',20,'NClones',6000);
% [NClones,NClones_ci95up,NClones_ci95dn] = calculate_CloneDens(timepoints,cloneSizes,showCI,sampling);

%% Default parameter values:
% Fill in unset optional values:
if (nargin < 3)
    showCI = 0;
    sampling = struct('NSubsets',20,'NClones',6000);
end

NClones = [];
NClones_ci95up = [];
NClones_ci95dn = [];

%% Calculate number of surviving clones (clone density):
for aja = 1:length(timepoints)
    cloneSizes_2plot{1,aja}(:,1) = cloneSizes(1,aja,:);
    NClones(aja) = size(find(cloneSizes_2plot{1,aja}~=0),1);
end

%% Calculate plausible intervals between simulation runs (by clonal subsampling):
if showCI == 1
    % Subsampling with replacement:
    SizeSubsets = sampling.NClones .* ones(1,length(timepoints)); % ~6000 clones at earliest time point (to ~fit experimental data)
    cloneSizes_2plot_sampled = subsampling_clones_overtime(sampling.NSubsets,SizeSubsets,cloneSizes_2plot,0,timepoints);
    % Calculate number of surviving clones from subsamples:
    for aja = 1:length(timepoints)
        for eje = 1:sampling.NSubsets
            NClones_samples(eje,aja) = size(find(cloneSizes_2plot_sampled{eje,1}{1,aja}~=0),1);
        end
    end
    % Calculate 95% plausible intervals on the number of surviving clones:
    NClones_ci95up = quantile(NClones_samples,0.975,1);
    NClones_ci95dn = quantile(NClones_samples,0.025,1);
end

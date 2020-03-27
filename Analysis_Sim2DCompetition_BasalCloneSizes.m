%% SCRIPT USED FOR ANALYZING SIMULATED MUTANT CLONE COMPETITION DYNAMICS UNDER DIFFERENT SCENARIOS:
% Mutant competition dynamics are simulated under a given experimental
% protocol and NCF-model parameter conditions, which can be chosen by the
% user.
% Properties of basal clone sizes and mutant cell populations are plotted.

%% SELECTION OF EXPERIMENTAL PROTOCOL / MODEL CONDITIONS TO BE SIMULATED:
% selectProtocol:
%   'DEN_monoMut'::  DEN induction & labelling at t=0 (mutants assigned SAME fitness: resilience to differentiate)
%   'DEN_polyMut'::  DEN induction & labelling at t=0 (mutants assigned DIFFERENT fitness: resilience to differentiate)
%   'CTL'::          Control conditions with no DEN induction, just labelling at t=0 (all cells are WT)
%   'DEN_Mml'::      DEN induction at time=0 (mutants assigned DIFFERENT fitness) followed by Mml* induction after a lag time
%   'Mml_Ind'::      Mml* induction at time=0 followed by a random induction (fluorescent label) after a lag time
%   'Custom'::       Set by user along with the parameters below in the code)

selectProtocol = 'DEN_monoMut'; % Protocol of choice from the ones above

%%%%%%%% ProtocolType 'Simple' %%%%%%%%
%   label.
%     |
%  -|DEN|-------------------
%    t=0

%%%%%% ProtocolType 'DEN_IndMml' %%%%%%
%           Mml*
%            |
%  -|DEN|-------------------
%    t=0    t=2m

%%%%%% ProtocolType 'IndMml_Ind' %%%%%%
%            label.
%              |
%  -|Mml*|-------------------
%    t=0      t=3m

%% SELECTION OF LATTICE & NCF-MODEL PARAMETER CONDITIONS:
% General lattice simulation parameters:
lattice.Dim = 200; % 2D grid size (number of cells per dimension of the square lattice)
lattice.Neigh = 6; % Neighborhood geometry (4, 6 or 8 neighbors per cell)
timelim = 61; % Simulation time span post-labelling (weeks) (DEN treatment is considered instantaneous)
nval = 488; % Number of (regularly spaced) time points for evaluation
indiv = 1; % Number of independent runs of tissue dynamics (usually fixed to 1)

% Choice of NCF-model parameter values:
if (strcmpi('DEN_monoMut',selectProtocol) || ...
    strcmpi('DEN_polyMut',selectProtocol) || ...
    strcmpi('CTL',selectProtocol) || ...
    strcmpi('DEN_IndMml',selectProtocol) || ...
    strcmpi('IndMml_Ind',selectProtocol))
    % Load default (preset) parameter values:
    ParamSet = SelectModelParamVal(selectProtocol);

elseif (strcmpi('Custom',selectProtocol))
    % Specific, customized parameter values:
    ParamSet.ProtocolType = 'Simple';
    ParamSet.modeMut = 'polyMut';
    ParamSet.Lambda = 0.5; %(/week)
    ParamSet.freqMut0 = 0.2;
    ParamSet.fitnessMutShape = 6;
    ParamSet.freqLabel = 0.05;
end

%% COMPUTATIONAL SIMULATION OF CLONE COMPETITION IN 2D UNDER GIVEN PARAMETERS
% Basal cells are simulated according to the single-progenitor (SP) model
% paradigm under the specific given conditions, and cellular population
% size properties retrieved.
% In this script, a CTL case is run along with the experimental mutagenesis
% condition (TTM = treatment) set above, for the purpose of clonal feature
% comparison between both scenarios.

% Selection of MonteCarlo simulator (depending on the complexity of protocol type):
if strcmpi('Simple',ParamSet.ProtocolType)
    [TTM_nx1_count,TTM_nx2_count,TTM_nx_basal,TTM_ntime,TTM_ALL_x_Type,TTM_ALL_x_Clone,TTM_ALL_x_Label] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics(timelim, ParamSet.Lambda, ParamSet.freqLabel, ParamSet, lattice, nval, indiv);
    ParamSet.freqMut0 = 0; % Run a CTL case (i.e. no DEN mutagenesis, freqMut0=0) under the same parameter conditions:
    [CTL_nx1_count,CTL_nx2_count,CTL_nx_basal,CTL_ntime,CTL_ALL_x_Type,CTL_ALL_x_Clone,CTL_ALL_x_Label] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics(timelim, ParamSet.Lambda, ParamSet.freqLabel, ParamSet, lattice, nval, indiv);
elseif strcmpi('DEN_IndMml',ParamSet.ProtocolType)
    [TTM_nx1_count,TTM_nx2_count,TTM_nx_basal,TTM_ntime,TTM_ALL_x_Type,TTM_ALL_x_Clone,TTM_ALL_x_Label,TTM_fitnessMut,TTM_nxp_basal] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics_challenge(timelim,ParamSet.Lambda,ParamSet.freqLabel,ParamSet.ProtocolType,ParamSet,lattice,nval,indiv);
    ParamSet.freqMut0 = 0; % Run a CTL case (i.e. no DEN mutagenesis, freqMut0=0, just Mml* induction) under the same parameter conditions:
    [CTL_nx1_count,CTL_nx2_count,CTL_nx_basal,CTL_ntime,CTL_ALL_x_Type,CTL_ALL_x_Clone,CTL_ALL_x_Label,CTL_fitnessMut,CTL_nxp_basal] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics_challenge(timelim,ParamSet.Lambda,ParamSet.freqLabel,ParamSet.ProtocolType,ParamSet,lattice,nval,indiv);
elseif strcmpi('IndMml_Ind',ParamSet.ProtocolType)
    [TTM_nx1_count,TTM_nx2_count,TTM_nx_basal,TTM_ntime,TTM_ALL_x_Type,TTM_ALL_x_Clone,TTM_ALL_x_Label,TTM_fitnessMut,TTM_nxp_basal] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics_challenge(timelim,ParamSet.Lambda,ParamSet.freqLabel,ParamSet.ProtocolType,ParamSet,lattice,nval,indiv);
    ParamSet.freqLabel = 0; % Run a CTL case (i.e. no Mml* induction, freqLabel=0, just subsequent labelling) under the same parameter conditions:
    [CTL_nx1_count,CTL_nx2_count,CTL_nx_basal,CTL_ntime,CTL_ALL_x_Type,CTL_ALL_x_Clone,CTL_ALL_x_Label,CTL_fitnessMut,CTL_nxp_basal] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics_challenge(timelim,ParamSet.Lambda,ParamSet.freqLabel,ParamSet.ProtocolType,ParamSet,lattice,nval,indiv);    
end

% Save data (workspace variables) into file 'myTest.mat' into ./Datasets folder:
if ~exist('Datasets', 'dir') % makes new directory if this does not exist in pwd
   mkdir('Datasets')
end
save ./Datasets/myTest.mat

%% PLOT DISPLAY SETTINGS:

% Select whether to calculate and plot plausible intervals on simulated estimates
% (plausible intervals reflect the level of uncertainty given a limited clone sampling, similar to the actual experimental one)
% (so far this option is only available - contrasted with realistic experimental sample sizes - in the case of 'Simple' protocols)
showCI = 0; % Show plausible intervals on model outcome ( 0=NO | 1=YES )  (be aware of the significant time consumption)

%% CALCULATE AND PLOT AVERAGE BASAL CLONE SIZE OVER TIME:

% Plausible interval settings:
if showCI == 1 % A No. of clones aprox. equivalent to the No. of experimental ones at time 0 is sampled for plausible interval calculations
    sampling = struct('NSubsets',20,'NClones',6000); % parameters of random permutation sampling for plausible interval calculation
else
    sampling = struct(); % void parameters (plausible interval not requested)
end

% Calculate:
[avgCloneSize_CTL.avg,avgCloneSize_CTL.ci95up,avgCloneSize_CTL.ci95dn] = calculate_AvgCloneSize(CTL_ntime,CTL_nx_basal,showCI,sampling);
[avgCloneSize_TTM.avg,avgCloneSize_TTM.ci95up,avgCloneSize_TTM.ci95dn] = calculate_AvgCloneSize(TTM_ntime,TTM_nx_basal,showCI,sampling);

% Plot (with or without plausible intervals):
plot_AvgCloneSize(CTL_ntime,avgCloneSize_CTL,showCI,struct('Color',[0 0.45 0.74],'Leyend','CTL'))
plot_AvgCloneSize(TTM_ntime,avgCloneSize_TTM,showCI,struct('Color',[0.85 0.33 0.01],'Leyend','Treat.'))

%% CALCULATE AND PLOT CLONE DENSITY OVER TIME:

% Plausible interval settings:
if showCI == 1 % A No. of clones aprox. equivalent to the No. of experimental ones at time 0 is sampled for plausible interval calculations
    sampling = struct('NSubsets',20,'NClones',6000); % parameters of random permutation sampling for plausible interval calculation
else
    sampling = struct(); % void parameters (plausible interval not requested)
end

% Calculate:
[CTL_NClones.avg,CTL_NClones.ci95up,CTL_NClones.ci95dn] = calculate_CloneDens(CTL_ntime,CTL_nx_basal,showCI,sampling);
[TTM_NClones.avg,TTM_NClones.ci95up,TTM_NClones.ci95dn] = calculate_CloneDens(TTM_ntime,TTM_nx_basal,showCI,sampling);

% Plot (with or without plausible intervals):
plot_CloneDens(CTL_ntime,CTL_NClones,showCI,struct('Color',[0 0.45 0.74],'Leyend','CTL'))
plot_CloneDens(TTM_ntime,TTM_NClones,showCI,struct('Color',[0.85 0.33 0.01],'Leyend','Treat.'))

%% CALCULATE AND PLOT BASAL CLONE SIZE DISTRIBUTIONS OVER TIME:
% Cumulative clone size distributions are retrieved, normalized by the
% average clone size at each time point.

% Plausible interval settings:
if showCI == 1 % A No. of clones equivalent to the No. of experimental ones at each time point is sampled for plausible interval calculations
    sampling_CTL = struct('NSubsets',20,'NClones',[11552 15865 4152 2474 3485]); % parameters of random permutation sampling in CTL condition
    sampling_TTM = struct('NSubsets',20,'NClones',[15092  5682  539  281  188]); % parameters of random permutation sampling in DEN condition
    time2calc = [10 1*30 3*30 6*30 12*30]./7; % time points to retrieve data from (weeks) must fit those from the experiments for plausible interval calculations
else
    sampling_CTL = struct(); % void parameters (plausible interval not requested)
    sampling_TTM = struct(); % void parameters (plausible interval not requested)
    time2calc = [10 1*30 3*30 6*30 12*30]./7; % time points to retrieve data from (weeks) can be changed freely.
end

% Calculate:
[CUM_nfreq_CTL.full,CUM_nfreq_CTL.centre,CUM_nfreq_CTL.ci95up,CUM_nfreq_CTL.ci95dn,avgCloneSize_CTL] = calculate_CloneSizeDist(CTL_ntime,time2calc,CTL_nx_basal,showCI,sampling_CTL);
[CUM_nfreq_TTM.full,CUM_nfreq_TTM.centre,CUM_nfreq_TTM.ci95up,CUM_nfreq_TTM.ci95dn,avgCloneSize_TTM] = calculate_CloneSizeDist(TTM_ntime,time2calc,TTM_nx_basal,showCI,sampling_TTM);

% Display settings:
time2plot = [10 1*30 3*30 6*30 12*30]./7; % time points for plotting (weeks) (always as many elements as time2calc or less)

% Plot (with or without plausible intervals):
plot_CloneSizeDist(time2calc,time2plot,CUM_nfreq_CTL,showCI,avgCloneSize_CTL,struct('Color',[0 0.45 0.74],'Leyend','CTL'))
plot_CloneSizeDist(time2calc,time2plot,CUM_nfreq_TTM,showCI,avgCloneSize_TTM,struct('Color',[0.85 0.33 0.01],'Leyend','Treat.'))

%% PLOT 2D VIEWS OF THE EPITHELIAL CLONES OVER TIME:

% Display settings:
time2plot = [10 1*30 3*30 6*30 12*30]./7; % time points to be plotted (weeks)
ZoomIn = 3; % zoom on the lattice view (e.g. 3x shows 1/9 of the total area)

% Plot:
plot_2Dview(CTL_ntime,time2plot,CTL_ALL_x_Clone,ZoomIn,struct('row',1,'rowSpan',2,'Leyend','CTL'))
plot_2Dview(TTM_ntime,time2plot,TTM_ALL_x_Clone,ZoomIn,struct('row',2,'rowSpan',2,'Leyend','Treat.'))

%% PLOT FRACTION OF MUTANT CELLS IN THE BASAL LAYER OVER TIME:

% Calculate: (nx2=all mutant cell types, which are then divided by all populations)
CTL_FreqMut = CTL_nx2_count ./ (CTL_nx1_count+CTL_nx2_count) .*100;
TTM_FreqMut = TTM_nx2_count ./ (TTM_nx1_count+TTM_nx2_count) .*100;

% Plot:
plot_FreqMut(CTL_ntime,CTL_FreqMut,struct('Color',[0 0.45 0.74],'Leyend','CTL'))
plot_FreqMut(TTM_ntime,TTM_FreqMut,struct('Color',[0.85 0.33 0.01],'Leyend','Treat.'))

%% CALCULATE AND PLOT FRACTION OF LABELLED CELLS OVER TIME:

% Calculate: (sum all labelled positions at a given time and divide by the full lattice dimensions)
CTL_FreqLabel = sum(sum(CTL_ALL_x_Label(1,:,:,:),4),3) ./ size(CTL_nx_basal,3) .*100;
TTM_FreqLabel = sum(sum(TTM_ALL_x_Label(1,:,:,:),4),3) ./ size(TTM_nx_basal,3) .*100;

% Plot:
plot_FreqLabel(CTL_ntime,CTL_FreqLabel,struct('Color',[0 0.45 0.74],'Leyend','CTL'))
plot_FreqLabel(TTM_ntime,TTM_FreqLabel,struct('Color',[0.85 0.33 0.01],'Leyend','Treat.'))

%% LAUNCHING SIMULATIONS UNDER THE 'Simple' PROTOCOL TYPE (monoMut assumption)
% Parameter values:
timelim = 52; % (weeks)
Lambda = 0.8; % replacement rate (/week)
freqLabel = 0.02; % initial frequency of labelled basal cells
ParamVal.ProtocolType = 'Simple';
ParamVal.modeMut = 'monoMut';
ParamVal.freqMut0 = 0.02; % frequency of mutant cells induced by mutagenesis (%)
ParamVal.fitnessMut = 0.5; % fitness gain for mutant cells induced by muatagenesis (all the same)
lattice.Dim = 100; % lattice size
lattice.Neigh = 6; % Neighborhood geometry (number of cell neighbors)
nval = 416; % No. of time points where data is collected
indiv = 1;

% Simulations:
[TTM_nx1_count,TTM_nx2_count,TTM_nx_basal,TTM_ntime,TTM_ALL_x_Type,TTM_ALL_x_Clone,TTM_ALL_x_Label,TTM_fitnessMut] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics(timelim,Lambda,freqLabel,ParamVal,lattice,nval,indiv);
ParamVal.freqMut0 = 0; % Run a CTL case (i.e. no DEN mutagenesis, freqMut0=0) under the same parameter conditions:
[CTL_nx1_count,CTL_nx2_count,CTL_nx_basal,CTL_ntime,CTL_ALL_x_Type,CTL_ALL_x_Clone,CTL_ALL_x_Label,CTL_fitnessMut] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics(timelim,Lambda,freqLabel,ParamVal,lattice,nval,indiv);


%% LAUNCHING SIMULATIONS UNDER THE 'Simple' PROTOCOL TYPE (polyMut assumption)
% Parameter values:
timelim = 52; % (weeks)
Lambda = 0.8; % replacement rate (/week)
freqLabel = 0.02; % initial frequency of labelled basal cells
ParamVal.ProtocolType = 'Simple';
ParamVal.modeMut = 'polyMut';
ParamVal.freqMut0 = 0.25; % frequency of mutant cells induced by mutagenesis (%)
ParamVal.fitnessMutShape = 6; % shape parameter of the distribution of fitness gain values for mutant cells induced by mutagenesis
lattice.Dim = 100; % lattice size
lattice.Neigh = 6; % Neighborhood geometry (number of cell neighbors)
nval = 416; % No. of time points where data is collected
indiv = 1;

% Simulations:
[TTM_nx1_count,TTM_nx2_count,TTM_nx_basal,TTM_ntime,TTM_ALL_x_Type,TTM_ALL_x_Clone,TTM_ALL_x_Label,TTM_fitnessMut] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics(timelim,Lambda,freqLabel,ParamVal,lattice,nval,indiv);
ParamVal.freqMut0 = 0; % Run a CTL case (i.e. no DEN mutagenesis, freqMut0=0) under the same parameter conditions:
[CTL_nx1_count,CTL_nx2_count,CTL_nx_basal,CTL_ntime,CTL_ALL_x_Type,CTL_ALL_x_Clone,CTL_ALL_x_Label,CTL_fitnessMut] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics(timelim,Lambda,freqLabel,ParamVal,lattice,nval,indiv);


%% LAUNCHING SIMULATIONS UNDER THE 'DEN_IndMml' PROTOCOL TYPE
% Parameter values:
timelim = 52+9; % to allow enough time (induction occurs after DEN)
Lambda = 0.8; % replacement rate (/week)
ParamVal.ProtocolType = 'DEN_IndMml';
ParamVal.freqLabel = 0.002; % frequency of Mml*-GFP induction (%)
ParamVal.freqMut0 = 0.25; % frequency of mutant cells induced by mutagenesis (%)
ParamVal.fitnessMutShape = 6; % shape parameter of the distribution of fitness gain values for mutant cells induced by mutagenesis
ParamVal.lagTime = 2/12*365/7; % time of Mml* induction post-DEN (weeks)
ParamVal.fitnessMut = 0.7; % fitness gain of Mml* clones
lattice.Dim = 100; % lattice size
lattice.Neigh = 6; % Neighborhood geometry (number of cell neighbors)
nval = 416; % No. of time points where data is collected
indiv = 1;

% Simulations:
[TTM_nx1_count,TTM_nx2_count,TTM_nx_basal,TTM_ntime,TTM_ALL_x_Type,TTM_ALL_x_Clone,TTM_ALL_x_Label,TTM_fitnessMut,TTM_nxp_basal] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics_challenge(timelim,Lambda,ParamVal.freqLabel,ParamVal.ProtocolType,ParamVal,lattice,nval,indiv);
ParamVal.freqMut0 = 0; % Run a CTL case (i.e. no DEN mutagenesis, freqMut0=0, just Mml* induction) under the same parameter conditions:
[CTL_nx1_count,CTL_nx2_count,CTL_nx_basal,CTL_ntime,CTL_ALL_x_Type,CTL_ALL_x_Clone,CTL_ALL_x_Label,CTL_fitnessMut,CTL_nxp_basal] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics_challenge(timelim,Lambda,ParamVal.freqLabel,ParamVal.ProtocolType,ParamVal,lattice,nval,indiv);


%% LAUNCHING SIMULATIONS UNDER THE 'IndMml_Ind' PROTOCOL TYPE
% Parameter values:
timelim = 52+9; % to allow enough time (second induction occurs after Mml* induction)
Lambda = 0.8; % replacement rate (/week)
ParamVal.ProtocolType = 'IndMml_Ind';
ParamVal.freqLabel = 0.002; % frequency of Mml*-GFP induction (%)
ParamVal.fitnessMut = 0.7; % fitness gain of Mml* clones
ParamVal.lagTime = 3/12*365/7; % time of Confetti induction post-Mml* induction (weeks)
ParamVal.freqLabel2 = 0.1; % frequency of Confetti-labelling induction (%)        
lattice.Dim = 100; % lattice size
lattice.Neigh = 6; % Neighborhood geometry (number of cell neighbors)
nval = 416; % No. of time points where data is collected
indiv = 1;

% Simulations:
[TTM_nx1_count,TTM_nx2_count,TTM_nx_basal,TTM_ntime,TTM_ALL_x_Type,TTM_ALL_x_Clone,TTM_ALL_x_Label,TTM_fitnessMut,TTM_nxp_basal] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics_challenge(timelim,Lambda,ParamVal.freqLabel,ParamVal.ProtocolType,ParamVal,lattice,nval,indiv);
ParamVal.freqLabel = 0; % Run a CTL case (i.e. no Mml* induction, freqLabel=0, just subsequent labelling) under the same parameter conditions:
[CTL_nx1_count,CTL_nx2_count,CTL_nx_basal,CTL_ntime,CTL_ALL_x_Type,CTL_ALL_x_Clone,CTL_ALL_x_Label,CTL_fitnessMut,CTL_nxp_basal] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics_challenge(timelim,Lambda,ParamVal.freqLabel,ParamVal.ProtocolType,ParamVal,lattice,nval,indiv);    


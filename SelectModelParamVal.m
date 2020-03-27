function [ParamVal] = SelectModelParamVal(selectProtocolTypeMode)
%% SELECT MODEL PARAMETER VALUES DEPENDING ON PROTOCOL OF INTEREST:
% A collection of some preset parameter values used in the manuscript to
% simulate different experimental protocol conditions is provided.
% The function helps selecting the specific set of interest.

% from Colom et al, 2020

%% Input:
% selectProtocolTypeMode: string; type of simulated experimental condition to load parameter values for

%% Output:
% ParamVal: structure containing the parameter values of interest for that simulated condition

%% Parameter selection:
switch selectProtocolTypeMode

    case 'DEN_monoMut'
        %% SUITABLE PARAMETER VALUES FOR: DEN (monoMut)
        ParamVal.ProtocolType = 'Simple';
        ParamVal.modeMut = 'monoMut';
        ParamVal.Lambda = 0.8; % replacement rate (/week)
        ParamVal.freqMut0 = 0.02; % frequency of mutant cells induced by mutagenesis (%)
        ParamVal.fitnessMut = 0.5; % fitness gain for mutant cells induced by muatagenesis (all the same)
        ParamVal.freqLabel = 0.02; % initial frequency of labelled basal cells

    case 'DEN_polyMut'
        %% SUITABLE PARAMETER VALUES FOR: DEN (polyMut)
        ParamVal.ProtocolType = 'Simple';
        ParamVal.modeMut = 'polyMut';
        ParamVal.Lambda = 0.8; % replacement rate (/week)
        ParamVal.freqMut0 = 0.25; % frequency of mutant cells induced by mutagenesis (%)
        ParamVal.fitnessMutShape = 6; % shape parameter of the distribution of fitness gain values for mutant cells induced by mutagenesis
        ParamVal.freqLabel = 0.02; % initial frequency of labelled basal cells

    case 'CTL'
        %% SUITABLE PARAMETER VALUES FOR: CTL
        ParamVal.ProtocolType = 'Simple';
        ParamVal.modeMut = 'monoMut';
        ParamVal.Lambda = 0.8; % replacement rate (/week)
        ParamVal.freqMut0 = 0; % no mutagenesis is induced
        ParamVal.fitnessMut = 0.5; % ignored parameter
        ParamVal.freqLabel = 0.02; % initial frequency of labelled basal cells
        
    case 'DEN_Mml'
        %% SUITABLE PARAMETER VALUES FOR: DEN + Mml*
        ParamVal.ProtocolType = 'DEN_IndMml';
        ParamVal.Lamda = 0.8; % replacement rate (/week)
        ParamVal.freqMut0 = 0.25; % frequency of mutant cells induced by mutagenesis (%)
        ParamVal.fitnessMutShape = 6; % shape parameter of the distribution of fitness gain values for mutant cells induced by mutagenesis
        ParamVal.lagTime = 2/12*365/7; % time of Mml* induction post-DEN (weeks)
        ParamVal.freqLabel = 0.002; % frequency of Mml*-GFP induction (%)
        ParamVal.fitnessMut = 0.7; % fitness gain of Mml* clones
        
    case 'Mml_Ind'
        %% SUITABLE PARAMETER VALUES FOR: Mml* + Confetti induction
        ParamVal.ProtocolType = 'IndMml_Ind';
        ParamVal.Lambda = 0.8; % replacement rate (/week)
        ParamVal.freqLabel = 0.002; % frequency of Mml*-GFP induction (%)
        ParamVal.fitnessMut = 0.7; % fitness gain of Mml* clones
        ParamVal.lagTime = 3/12*365/7; % time of Confetti induction post-Mml* induction (weeks)
        ParamVal.freqLabel2 = 0.1; % frequency of Confetti-labelling induction (%)

end

function [nx1_count,nx2_count,nx_basal,ntime,ALL_x_Type,ALL_x_Clone,ALL_x_Label,fitnessMut,nxp_basal] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics_challenge(timelim,Lambda,freqLabel,ProtocolType,ProtocolVars,lattice,nval,indiv)
%% Lattice-based Monte Carlo simulator of cell competition dynamics in poly-mutated tissue (version comprising complex protocols)
% Competing epithelial progenitor cells are simulated in a 2-dimensional
% grid representing the basal proliferative compartment. These cells divide
% at random replacing a neighbor. Wild-type and mutant progenitors are
% considered, with different propensity to differentiate (be replaced).
% Individual clones are tracked over the time course.

% This version is particularly suited for simulating complex, two-stage
% protocols involving initial mutagenesis and a late induction (intervention)
% in the middle of the period of chase. This is the case of a protocol
% ('DEN_IndMml') where random mutagenesis is assumed to occur at time 0 and
% a particular mutation is induced later on after a lag time; or a protocol
% ('IndMml_Ind') where a particular mutation is induced at time 0 and
% random labelling occurs later on after a lag time. These scenarios
% simulate the experimental conditions where Mml*-GFP is induced after a
% period of DEN-driven mutagenesis, and the experiments where Confetti-
% lineage tracing is performed in mice upon Mml*-GFP induction, respectively.

% from Colom et al, 2020

%% Input:
% REQUIRED:
% timelim: horizontal vector containing the desired time span (expressed in weeks)
% Lambda: average replacement rate (/week)
% freqLabel: frequency of labelling-induced basal cells (whether RYFP-, Confetti- or Mml*-GFP cells)
% ProtocolType: string indicating the protocol type
    % 'Simple': Mutagenesis at time 0, before chase
    % 'IndMml_Ind': Mml*-GFP induction at time 0; random labelling during chase
    % 'DEN_IndMml': Mutagenesis at time 0; Mml*-GFP induction during chase
% ProtocolVars: structure containing protocol-related parameters (depend on ProtocolType):
    %%%%%%%% 'Simple' protocol %%%%%%%%:
    % ProtocolVars == struct{modeMut==monoMut, freqMut0, fitnessMut} -or- struct{modeMut==polyMut, freqMut0, fitnessMutShape}
        % modeMut: compaction mode (either 'monoMut' or 'polyMut', depending on whether all mutant cells are given same fitness gain or not)
        % freqMut0: initial frequency of mutant cells (value between 0 and 1)
        % fitnessMutShape: Shape parameter of the Gamma distribution for mutant fitness gain values. Only in 'polyMut' mode.
        % fitnessMut: Single fitness gain value assigned to all mutants. Only in 'monoMut' mode.
    %%%%%%%% 'IndMml_Ind' protocol %%%%%%%%:
    % ProtocolVars == struct{fitnessMut, lagTime, freqLabel2}
        % fitnessMut: Single fitness gain value assigned to Mml* mutant genotype.
        % lagTime: time when random labelling occurs after Mml*-GFP induction
        % freqLabel2: frequency of random labelling in this second stage
    %%%%%%%% 'DEN_IndMml' protocol %%%%%%%%:
    % ProtocolVars == struct{freqMut0, fitnessMutShape, lagTime, fitnessMut}
        % freqMut0: initial frequency of mutant cells caused by mutagenesis (value between 0 and 1)
        % fitnessMutShape: Shape parameter of the Gamma distribution for fitness gain values of mutants induced by mutagenesis.
        % lagTime: time when Mml*-GFP induction occurs after mutagenesis
        % fitnessMut: Single fitness gain value assigned to Mml* mutant genotype.

% OPTIONAL:
% lattice: structure containing lattice-related parameters:
    % struct{Dim==100, Neigh}    
        % Dim: number of cells per dimension of the square lattice, e.g. 100x100
        % Neigh: number of neighbors for a given basal cell (geometry) (value = 4, 6 or 8)
% nval: number of (regularly spaced) time points for evaluation
% indiv: number of independent runs of tissue dynamics

%% Output:
% nx1_count: total no. of wild-type cells over time
% nx2_count: total no. of mutant cells over time
% nx_basal: mxn matrix containing sizes (no. of basal cells) of all m initial clones over the n time points
% ntime: horizontal vector of the n time points collected
% ALL_x_Type: matrix reporting the genotype for each cell location at each time point
% ALL_x_Clone: matrix reporting the clone ID for each cell location at each time point
% ALL_x_Label: matrix reporting whether each cell location is labelled or not at each time point
% fitnessMut: px1 column vector containing the fitness gain values for the p different genotypes (1st element is WT)
% nxp_basal: mxn matrix containing sizes (no. of basal cells) of all subclones induced at 2nd stage over the n time points (case of 'IndMml_Ind' or 'DEN_IndMml' protocol)

%% Example:
% timelim = 52+9; %(weeks) e.g. to comprise at least 1y after induction (which may occur after mutegenesis)
% Lambda = 0.8; %(/week)
% freqLabel = 0.002;
% ProtocolType = 'DEN_IndMml';
% ProtocolVars.freqMut0 = 0.25;
% ProtocolVars.fitnessMutShape = 6;
% ProtocolVars.lagTime = 2/12*365/7 % (weeks)
% ProtocolVars.fitnessMut = 0.7;
% lattice.Dim = 100;
% lattice.Neigh = 6;
% nval = 488;
% indiv = 1;
% [nx1_count,nx2_count,nx_basal,ntime,ALL_x_Type,ALL_x_Clone,ALL_x_Label,fitnessMut,nxp_basal] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics_challenge(timelim,Lambda,freqLabel,ProtocolType,ProtocolVars,lattice,nval,indiv);

%% Default parameter values:
tic
% Check number of inputs:
if nargin < 5
    error('MonteCarloSim:TooFewInputs', ...
        'requires at least 5 inputs');
end
% Fill in unset optional values:
if nargin < 8
    lattice = struct('Dim',100,'Neigh',6);
    nval=416; % no. of registered time points (for variable recording)
    indiv = 1; % no. of 2D-tissue samples
end

% Initial definition of parameters:
parcialtime=timelim./nval;
ntime = ones(indiv,1) * [0:nval].*parcialtime;
ALL_x_Type = zeros(indiv,nval+1,lattice.Dim,lattice.Dim);
ALL_x_Clone = zeros(indiv,nval+1,lattice.Dim,lattice.Dim);
ALL_x_subClone = zeros(indiv,nval+1,lattice.Dim,lattice.Dim);
ALL_x_Label = zeros(indiv,nval+1,lattice.Dim,lattice.Dim);
nx1_count = zeros(indiv,nval+1);
nx2_count = zeros(indiv,nval+1);
nx_basal = zeros(indiv,nval+1,lattice.Dim*lattice.Dim);
nxp_basal = zeros(indiv,nval+1,lattice.Dim*lattice.Dim);

% Vector of neighbour locations respect to any particular target cell (neighbourhood geometry):
if lattice.Neigh == 4
    % Square geometry, 4 neighbors (Von Neumann's layout):
    neigh_pos = [-1 0 0 +1;... % row displacement
        0 -1 +1 0]; % col displacement
elseif lattice.Neigh == 8
    % Square geometry, 8 neighbors (Moore layout):
    neigh_pos = [-1 -1 -1 0 0 +1 +1 +1;... % row displacement
        -1 0 +1 -1 +1 -1 0 +1]; % col displacement
else
    % Hexagonal geometry, 6 neighbors:
    neigh_pos = [-1 -1 0 0 +1 +1;... % row displacement
        0 +1 -1 +1 -1 0]; % col displacement    
end

%% ITERATION FOR DIFFERENT INDEPENDENT RUNS OF TISSUE DYNAMICS
for it=1:indiv

    it;

    % Define initial cell attributes (depending on the protocol of choice):
    x_Clone = reshape([1:lattice.Dim*lattice.Dim],lattice.Dim,lattice.Dim);% each initial basal cell defines a single clone
    switch ProtocolType
        % depending on the protocol, both mutagenesis & labelling (Simple & IndMml_Ind) or just mutagenesis (DEN_IndMml) occur at time 0
        case 'Simple'
            % induce random mutagenesis:
            x_Type = double( rand(lattice.Dim,lattice.Dim) < (ProtocolVars.freqMut0) ); % 0 means position occupied by WT cell | 1 means position occupied by Mut cell.
            if ProtocolVars.modeMut == 'polyMut'
                x_Type(find(x_Type(:)>0)) = [1:size(find(x_Type(:)>0),1)]; % assign a different code to each individual mutant, i.e. 1,2,...,m
                fitnessMut = 1-gamrnd(ProtocolVars.fitnessMutShape,1/ProtocolVars.fitnessMutShape,[size(find(x_Type(:)>0),1) 1]); % assign random fitness values to each individual mutant (from an underlying Gamma dist)
            else %'monoMut'
                fitnessMut = ProtocolVars.fitnessMut; % assign a single fitness value to all mutants (x_Type already defined)
            end
            % induce random labelling:
            x_Label = rand(lattice.Dim,lattice.Dim) < (freqLabel); % 0 means non-labelled clone | 1 means labelled clone
        
        case {'IndMml_Ind'}
            % induce Mml* mutation & labelling:
            x_Type = double( rand(lattice.Dim,lattice.Dim) < (freqLabel) ); % 0 means position occupied by WT cell | 1 means position occupied by Mut cell.
            if max(x_Type(:))>0
                fitnessMut = ProtocolVars.fitnessMut; % assign a single fitness value to Mml mutants (x_Type already defined)
            else
                fitnessMut = []; % leave vector of fitness values empty
            end
            x_Label = x_Type; % OR = logical(x_Type) | label induction of Mml*
            x_subClone = zeros(lattice.Dim,lattice.Dim); % subclonal labelling induction occurs later
        
        case {'DEN_IndMml'}
            % induce random mutagenesis:
            x_Type = double( rand(lattice.Dim,lattice.Dim) < (ProtocolVars.freqMut0) ); % 0 means position occupied by WT cell | 1 means position occupied by Mut cell.
            x_Type(find(x_Type(:)>0)) = [1:size(find(x_Type(:)>0),1)]; % assign a different code to each individual mutant
            fitnessMut = 1-gamrnd(ProtocolVars.fitnessMutShape,1/ProtocolVars.fitnessMutShape,[size(find(x_Type(:)>0),1) 1]); % assign random fitness values to each individual mutant (from an underlying Gamma dist)
            % no labelling induction at time 0:
            x_Label = zeros(lattice.Dim,lattice.Dim); % label induction occurs later
            x_subClone = zeros(lattice.Dim,lattice.Dim); % label induction occurs later
    end

    % Initial variables:
    time=0;
    multiplo=0;
    ProtocolStage = 0;

    % Calculation of number of cells of each type at time 0:
    x1_count = size(find(x_Type==0),1); % total no. of WT cells
    x2_count = size(find(x_Type>=1),1); % total no. of Mut cells (no matter the particular mutation)

    % Save the locations and populations of cells at time 0:
    ALL_x_Type(it,multiplo+1,:,:) = x_Type(:,:);
    ALL_x_Clone(it,multiplo+1,:,:) = x_Clone(:,:);
    ALL_x_Label(it,multiplo+1,:,:) = x_Label(:,:);
    nx1_count(it,multiplo+1) = x1_count;
    nx2_count(it,multiplo+1) = x2_count;
    nx_basal(it,multiplo+1,:) = histc(x_Clone(:),[1:1:lattice.Dim*lattice.Dim]);
    ntime(it,multiplo+1) = time;
    multiplo = multiplo+1;

    % ITERATION FOR EACH SINGLE TISSUE:
    while (time <= timelim)

        % Second part of protocol (if any):
        try % (lagTime is not a variable in the 'Simple' protocol)
            if (time > ProtocolVars.lagTime) && (ProtocolStage == 0) && (strcmp(ProtocolType,'IndMml_Ind'))
                % induce random labelling:
                x_LabelUpdate = double(  rand(lattice.Dim,lattice.Dim) < (ProtocolVars.freqLabel2) ); % 0 means position occupied by original label type | 1 means position occupied by new label.
                x_LabelUpdated = double(x_Label) + (x_LabelUpdate.*2); % 0=WT | 1=Mml | 2=WT_labelled | 3=Mml_labelled
                % assign a different code to each individual new labelled clone within WT background (2)
                try x_Label(find(x_LabelUpdated(:)==2)) = [2:size(find(x_LabelUpdated(:)==2),1)+1]; end
                % assign a different code to each individual new labelled clone sitting on a Mml background (3)
                try x_Label(find(x_LabelUpdated(:)==3)) = [ size(find(x_LabelUpdated(:)==2),1)+2 : size(find(x_LabelUpdated(:)==2),1)+1 + size(find(x_LabelUpdated(:)==3),1) ]; end
                % assign a different code (subclone ID) to each individual new labelled cell:
                x_subClone(find(x_LabelUpdate(:)>0)) = [1:size(find(x_LabelUpdate(:)>0),1)];
                % set the status of second part of protocol as completed:
                ProtocolStage = 1;
                
            elseif (time > ProtocolVars.lagTime) && (ProtocolStage == 0) && (strcmp(ProtocolType,'DEN_IndMml'))
                % induce Mml* mutation & labelling:
                x_TypeUpdate = double( rand(lattice.Dim,lattice.Dim) < (freqLabel) ); % 0 means position occupied by original cell type | 1 means position occupied by new Mut cell.
                x_TypeUpdated = 2.*x_Type + x_TypeUpdate; % 0=WT | 1=Mml | 2,4,6,..=Mut | 3,5,7,..=Mut+Mml
                % assign a single code and single fitness value to all Mml mutants occurring on WT background (1)
                try x_Type(find(x_TypeUpdated(:)==1)) = size(fitnessMut,1)+1; end
                if size(find(x_TypeUpdated(:)==1),1)>0
                    fitnessMut = [fitnessMut; ProtocolVars.fitnessMut];
                end
                % assign a different code and fitness value to each individual Mml sitting on a previous mutant (3,5,7,..)
                try x_Type(find(rem(x_TypeUpdated(:),2)>0 & x_TypeUpdated(:)>1)) = [ size(fitnessMut,1)+1 : size(fitnessMut,1) + size(find(rem(x_TypeUpdated(:),2)>0 & x_TypeUpdated(:)>1),1) ]; end
                %fitnessMut = [fitnessMut; 0.5-gamrnd(20,0.3/20,[size(find(rem(x_TypeUpdated(:),2)>0 & x_TypeUpdated(:)>1),1) 1])]; % random fitness values to each individual Mml sitting on a previous mutant background (from an underlying Gamma dist)
                fitnessMut = [fitnessMut; ProtocolVars.fitnessMut.*ones(size(find(rem(x_TypeUpdated(:),2)>0 & x_TypeUpdated(:)>1),1),1)]; % same fitness value to all individual Mml hitting previous mutant clones
                x_Label = logical(x_TypeUpdate);
                % assign a different code (subclone ID) to each individual new labelled cell:
                x_subClone(find(x_Label(:)>0)) = [1:size(find(x_Label(:)>0),1)];
                % set the status of second part of protocol as completed:
                ProtocolStage = 1;
            end
        end

        % Random number generator:
        r1=rand; % for time progression
        r2=rand; % for selection of cell to divide
        r3=rand; % for selection of neighbour to differentiate (displacement)

        % Calculation of time for next event: (just division)
        pt = Lambda*lattice.Dim*lattice.Dim; % WT -> WT + WT | Mut -> Mut + Mut
        tau=-(1./pt)*log(r1);

        % Selection of individual target cell-to-divide:
        target_element = ceil(r2*lattice.Dim*lattice.Dim);
        [target_locY,target_locX] = ind2sub(size(x_Type),target_element);

        % Define neighbour locations:
        neigh_locY = target_locY + neigh_pos(1,:);
        neigh_locX = target_locX + neigh_pos(2,:);
        % correct neighbour positions falling out of the lattice boundaries... (periodic boundary counditions)
        neigh_locY(find(neigh_locY<1)) = lattice.Dim; neigh_locY(find(neigh_locY>lattice.Dim)) = 1;
        neigh_locX(find(neigh_locX<1)) = lattice.Dim; neigh_locX(find(neigh_locX>lattice.Dim)) = 1;

        % Selection of neighbour cell-to-differentiate:
        pDiff = zeros(length(neigh_pos),1);
        for aja = 1:length(neigh_pos)
            if x_Type(neigh_locY(aja),neigh_locX(aja)) == 0; pDiff(aja) = 1; %WT prob. to differentiate = 1 (default)
            else; pDiff(aja) = 1-fitnessMut( x_Type(neigh_locY(aja),neigh_locX(aja)) ,1); %Mut prob. to differentiate = 1-fitnessMut
            end
        end
        ptDiff = sum(pDiff);
        whoDiff = find(cumsum(pDiff)>(r3*ptDiff),1);

        % Apply changes (the divided cell occupies the space of the differentiated cell):
        x_Type(neigh_locY(whoDiff),neigh_locX(whoDiff)) = x_Type(target_locY,target_locX);
        x_Clone(neigh_locY(whoDiff),neigh_locX(whoDiff)) = x_Clone(target_locY,target_locX);
        x_Label(neigh_locY(whoDiff),neigh_locX(whoDiff)) = x_Label(target_locY,target_locX);
        if ( strcmp(ProtocolType,'IndMml_Ind') || strcmp(ProtocolType,'DEN_IndMml') )
            x_subClone(neigh_locY(whoDiff),neigh_locX(whoDiff)) = x_subClone(target_locY,target_locX);
        end

        % Update the populations of each cell type:
        x1_count = size(find(x_Type==0),1); % total no. of WT cells
        x2_count = size(find(x_Type>=1),1); % total no. of Mut cells

        % Calculate time to that event:
        time=time+tau;

        % Save the populations/properties of cells at certain time points:
        if ( (time >= multiplo*parcialtime) && (time ~= Inf) )
            if (time < (multiplo+1)*parcialtime)
                ALL_x_Type(it,multiplo+1,:,:) = x_Type(:,:);
                ALL_x_Clone(it,multiplo+1,:,:) = x_Clone(:,:);
                ALL_x_subClone(it,multiplo+1,:,:) = x_subClone(:,:);
                ALL_x_Label(it,multiplo+1,:,:) = x_Label(:,:);
                nx1_count(it,multiplo+1) = x1_count; nx2_count(it,multiplo+1) = x2_count;
                nx_basal(it,multiplo+1,:) = histc(x_Clone(:),[1:1:lattice.Dim*lattice.Dim]);
                nxp_basal(it,multiplo+1,:) = histc(x_subClone(:),[1:1:lattice.Dim*lattice.Dim]);
                ntime(it,multiplo+1) = time;
                multiplo=multiplo+1;
            else
                multiplo_pre = multiplo;
                while ( (time >= (multiplo+1)*parcialtime) && ((multiplo+1)*parcialtime < (timelim+parcialtime)) )
                    ALL_x_Type(it,multiplo+1,:,:) = ALL_x_Type(it,multiplo_pre,:,:);
                    ALL_x_Clone(it,multiplo+1,:,:) = ALL_x_Clone(it,multiplo_pre,:,:);
                    ALL_x_subClone(it,multiplo+1,:,:) = ALL_x_subClone(it,multiplo_pre,:,:);
                    ALL_x_Label(it,multiplo+1,:,:) = ALL_x_Label(it,multiplo_pre,:,:);
                    nx1_count(it,multiplo+1) = nx1_count(it,multiplo_pre); nx2_count(it,multiplo+1) = nx2_count(it,multiplo_pre);
                    past_Clone(:,:) = ALL_x_Clone(it,multiplo_pre,:,:);
                    past_subClone(:,:) = ALL_x_subClone(it,multiplo_pre,:,:);
                    nx_basal(it,multiplo+1,:) = histc(past_Clone(:),[1:1:lattice.Dim*lattice.Dim]);
                    nxp_basal(it,multiplo+1,:) = histc(past_subClone(:),[1:1:lattice.Dim*lattice.Dim]);
                    ntime(it,multiplo+1)= multiplo*parcialtime;
                    multiplo=multiplo+1;
                end
                if (time >= timelim+parcialtime)
                    ALL_x_Type(it,multiplo+1,:,:) = ALL_x_Type(it,multiplo,:,:);
                    ALL_x_Clone(it,multiplo+1,:,:) = ALL_x_Clone(it,multiplo,:,:);
                    ALL_x_subClone(it,multiplo+1,:,:) = ALL_x_subClone(it,multiplo,:,:);
                    ALL_x_Label(it,multiplo+1,:,:) = ALL_x_Label(it,multiplo,:,:);
                    nx1_count(it,multiplo+1) = nx1_count(it,multiplo); nx2_count(it,multiplo+1) = nx2_count(it,multiplo);
                    pre_Clone(:,:) = ALL_x_Clone(it,multiplo,:,:);
                    pre_subClone(:,:) = ALL_x_subClone(it,multiplo,:,:);
                    nx_basal(it,multiplo+1,:) = histc(pre_Clone(:),[1:1:lattice.Dim*lattice.Dim]);
                    nxp_basal(it,multiplo+1,:) = histc(pre_subClone(:),[1:1:lattice.Dim*lattice.Dim]);
                    ntime(it,multiplo+1) = multiplo*parcialtime;
                else
                    ALL_x_Type(it,multiplo+1,:,:) = x_Type(:,:);
                    ALL_x_Clone(it,multiplo+1,:,:) = x_Clone(:,:);
                    ALL_x_subClone(it,multiplo+1,:,:) = x_subClone(:,:);
                    ALL_x_Label(it,multiplo+1,:,:) = x_Label(:,:);
                    nx1_count(it,multiplo+1) = x1_count; nx2_count(it,multiplo+1) = x2_count;
                    nx_basal(it,multiplo+1,:) = histc(x_Clone(:),[1:1:lattice.Dim*lattice.Dim]);
                    nxp_basal(it,multiplo+1,:) = histc(x_subClone(:),[1:1:lattice.Dim*lattice.Dim]);
                    ntime(it,multiplo+1)=time;
                end
                multiplo=multiplo+1;
            end
        end

    end

end

toc

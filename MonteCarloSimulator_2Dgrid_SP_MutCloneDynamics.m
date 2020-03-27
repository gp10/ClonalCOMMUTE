function [nx1_count,nx2_count,nx_basal,ntime,ALL_x_Type,ALL_x_Clone,ALL_x_Label,fitnessMut] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics(timelim,Lambda,freqLabel,ProtocolVars,lattice,nval,indiv)
%% Lattice-based Monte Carlo simulator of cell competition dynamics in poly-mutated tissue
% Competing epithelial progenitor cells are simulated in a 2-dimensional
% grid representing the basal proliferative compartment. These cells divide
% at random replacing a neighbor. Wild-type and mutant progenitors are
% considered, with different propensity to differentiate (be replaced).
% Individual clones are tracked over the time course.

% This version considers a simple protocol where random mutagenesis and
% independent, random labelling are assumed to occur at time 0 before the
% period of chase, a scenario simulating experimental RYFP- or Confetti-
% lineage tracing in mice following treatment with mutagen DEN.

% from Colom et al, 2020

%% Input:
% REQUIRED:
% timelim: horizontal vector containing the desired time span (expressed in weeks)
% Lambda: average replacement rate (/week)
% freqLabel: initial frequency of labelled basal cells
% ProtocolVars: structure containing mutant population parameters:
    % struct{modeMut==monoMut, freqMut0, fitnessMut} -or- struct{modeMut==polyMut, freqMut0, fitnessMutShape}
        % modeMut: compaction mode (either 'monoMut' or 'polyMut', depending on whether all mutant cells are given same fitness gain or not)
        % freqMut0: initial frequency of mutant cells (value between 0 and 1)
        % fitnessMutShape: Shape parameter of the Gamma distribution for mutant fitness gain values. Only in 'polyMut' mode.
        % fitnessMut: Single fitness gain value assigned to all mutants. Only in 'monoMut' mode.

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

%% Example:
% timelim = 52; %(weeks)
% Lambda = 0.8; %(/week)
% freqLabel = 0.02;
% ProtocolVars.modeMut = 'monoMut';
% ProtocolVars.freqMut0 = 0.02;
% ProtocolVars.fitnessMut = 0.5;
% lattice.Dim = 300;
% lattice.Neigh = 6;
% nval = 416;
% indiv = 1;
% [nx1_count,nx2_count,nx_basal,ntime,ALL_x_Type,ALL_x_Clone,ALL_x_Label,fitnessMut] = MonteCarloSimulator_2Dgrid_SP_MutCloneDynamics(timelim,Lambda,freqLabel,ProtocolVars,lattice,nval,indiv);

%% Default parameter values:
tic
% Check number of inputs:
if nargin < 4
    error('MonteCarloSim:TooFewInputs', ...
        'requires at least 4 inputs');
end
% Fill in unset optional values:
if nargin < 7
    lattice = struct('Dim',100,'Neigh',6);
    nval=416; % no. of registered time points (for variable recording)
    indiv = 1; % no. of 2D-tissue samples
end

% Initial definition of parameters:
parcialtime=timelim./nval;
ntime = ones(indiv,1) * [0:nval].*parcialtime;
ALL_x_Type = zeros(indiv,nval+1,lattice.Dim,lattice.Dim);
ALL_x_Clone = zeros(indiv,nval+1,lattice.Dim,lattice.Dim);
ALL_x_Label = zeros(indiv,nval+1,lattice.Dim,lattice.Dim);
nx1_count = zeros(indiv,nval+1);
nx2_count = zeros(indiv,nval+1);
nx_basal = zeros(indiv,nval+1,lattice.Dim*lattice.Dim);

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

    % Define initial cell types (random type, depending on freqMut0, and random fitnessMut)
    x_Clone = reshape([1:lattice.Dim*lattice.Dim],lattice.Dim,lattice.Dim);% each initial basal cell defines a single clone
    x_Type = double( rand(lattice.Dim,lattice.Dim) < (ProtocolVars.freqMut0) ); % 0 means position occupied by WT cell | 1 means position occupied by Mut cell.
    if ProtocolVars.modeMut == 'polyMut'
        x_Type(find(x_Type(:)>0)) = [1:size(find(x_Type(:)>0),1)]; % assign a different code to each individual mutant
        fitnessMut = 1-gamrnd(ProtocolVars.fitnessMutShape,1/ProtocolVars.fitnessMutShape,[size(find(x_Type(:)>0),1) 1]); % assign random fitness values to each individual mutant (from an underlying Gamma dist)
    else %'monoMut'
        fitnessMut = ProtocolVars.fitnessMut; % assign a single fitness value to all mutants (x_Type already defined)
    end
    x_Label = rand(lattice.Dim,lattice.Dim) < (freqLabel); % 0 means non-labelled clone | 1 means labelled clone

    % Initial variables:
    time=0;
    multiplo=0;

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

        % Random number generator:
        r1=rand; % for time progression
        r2=rand; % for selection of cell to divide
        r3=rand; % for selection of neighbour cell to differentiate (displacement)

        % Calculation of time for next event: (just division)
        pt = Lambda*lattice.Dim*lattice.Dim; % WT -> WT + WT | Mut -> Mut + Mut
        tau=-(1./pt)*log(r1);

        % Selection of individual target cell-to-divide:
        target_element = ceil(r2*lattice.Dim*lattice.Dim);
        [target_locY,target_locX] = ind2sub(size(x_Type),target_element);

        % Define neighbour locations:
        neigh_locY = target_locY + neigh_pos(1,:);
        neigh_locX = target_locX + neigh_pos(2,:);
        % correct neighbour positions falling out of the lattice boundaries... (periodic boundary conditions)
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
                ALL_x_Label(it,multiplo+1,:,:) = x_Label(:,:);
                nx1_count(it,multiplo+1) = x1_count; nx2_count(it,multiplo+1) = x2_count;
                nx_basal(it,multiplo+1,:) = histc(x_Clone(:),[1:1:lattice.Dim*lattice.Dim]);
                ntime(it,multiplo+1) = time;
                multiplo=multiplo+1;
            else
                multiplo_pre = multiplo;
                while ( (time >= (multiplo+1)*parcialtime) && ((multiplo+1)*parcialtime < (timelim+parcialtime)) )
                    ALL_x_Type(it,multiplo+1,:,:) = ALL_x_Type(it,multiplo_pre,:,:);
                    ALL_x_Clone(it,multiplo+1,:,:) = ALL_x_Clone(it,multiplo_pre,:,:);
                    ALL_x_Label(it,multiplo+1,:,:) = ALL_x_Label(it,multiplo_pre,:,:);
                    nx1_count(it,multiplo+1) = nx1_count(it,multiplo_pre); nx2_count(it,multiplo+1) = nx2_count(it,multiplo_pre);
                    past_Clone(:,:) = ALL_x_Clone(it,multiplo_pre,:,:);
                    nx_basal(it,multiplo+1,:) = histc(past_Clone(:),[1:1:lattice.Dim*lattice.Dim]);
                    ntime(it,multiplo+1)= multiplo*parcialtime;
                    multiplo=multiplo+1;
                end
                if (time >= timelim+parcialtime)
                    ALL_x_Type(it,multiplo+1,:,:) = ALL_x_Type(it,multiplo,:,:);
                    ALL_x_Clone(it,multiplo+1,:,:) = ALL_x_Clone(it,multiplo,:,:);
                    ALL_x_Label(it,multiplo+1,:,:) = ALL_x_Label(it,multiplo,:,:);
                    nx1_count(it,multiplo+1) = nx1_count(it,multiplo); nx2_count(it,multiplo+1) = nx2_count(it,multiplo);
                    pre_Clone(:,:) = ALL_x_Clone(it,multiplo,:,:);
                    nx_basal(it,multiplo+1,:) = histc(pre_Clone(:),[1:1:lattice.Dim*lattice.Dim]);
                    ntime(it,multiplo+1) = multiplo*parcialtime;
                else
                    ALL_x_Type(it,multiplo+1,:,:) = x_Type(:,:);
                    ALL_x_Clone(it,multiplo+1,:,:) = x_Clone(:,:);
                    ALL_x_Label(it,multiplo+1,:,:) = x_Label(:,:);
                    nx1_count(it,multiplo+1) = x1_count; nx2_count(it,multiplo+1) = x2_count;
                    nx_basal(it,multiplo+1,:) = histc(x_Clone(:),[1:1:lattice.Dim*lattice.Dim]);
                    ntime(it,multiplo+1)=time;
                end
                multiplo=multiplo+1;
            end
        end

    end

end

toc

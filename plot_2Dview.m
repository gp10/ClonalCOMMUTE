function plot_2Dview(timepoints,time2plot,Clone_2Dmaps,ZoomIn,FigProp)
%% Plots spatial 2D views of basal-layer clones over time
% It plots a spatial (top-down) view of clones in the basal layer at time
% points of interest, i.e. clones are displayed in the 2D lattice color-coded
% according to initial clone IDs.

% from Colom et al, 2020

%% Input:
% timepoints: vector of time points in the original data (expressed in weeks)
% time2plot: vector of time points of interest to plot data from (expressed in weeks)
% Clone_2Dmaps: matrix containing the spatial maps of clone IDs at different time points (as retrieved from the lattice-based simulator)
    % format: mxnxoxp (m=No.repeats; n=No.time points; o=y-coordinate; p=x-coordinate) 
% ZoomIn: visual zoom on the grid (1 = full grid-size view)
% FigProp: structure containing output display properties
    % struct{rowSpan==2, row==1, Leyend=='CTL'}
        % rowSpan: total number of rows in the multi-panel figure (e.g. designed to span for the total number of conditions to be compared)
        % row: specific row number in the multi-panel figure where time-series data is to be plotted
        % Leyend: legend entry

%% Example:
% timepoints = [1:10];
% time2plot = [2,6];
% Clone_2Dmaps = rand(1,10,40,40).*2000;
% ZoomIn = 1;
% FigProp.rowSpan = 2; FigProp.row = 1; FigProp.Leyend = 'CTL';
% plot_2Dview(timepoints,time2plot,Clone_2Dmaps,ZoomIn,FigProp);

%% Figure panel selection:
% Make handle to a new figure if this does not exist yet:
if isempty(findobj( 'Type', 'Figure', 'Name', '2Dview'))
    % makes new figure (with white background color):
    figure('Color',[1 1 1])
    set(gcf,'Name','2Dview');
    setHeader = 1; % sets this is first data in the figure
else
    % retrieves existing one:
    figure(findobj( 'Type', 'Figure', 'Name', '2Dview'))
    setHeader = 0; % sets there were preceding data in the figure
end

%% Plot spatial 2D view of basal-layer clones over time:
% Settings before plotting:
latticeDim = size(Clone_2Dmaps,3);
myColClones = rand(latticeDim^2,3); % ensures each clone ID will be mapped to a different color

% Plot data in new row:
for aja = 1:length(time2plot)
    myTime = find(timepoints >= time2plot(aja),1); % position in the time vector
    Clone_2Dmaps_2plot(:,:) = Clone_2Dmaps(1,myTime, 1:round(latticeDim/ZoomIn), 1:round(latticeDim/ZoomIn)); % time-specific map of clone locations
    
    subplot(FigProp.rowSpan,length(time2plot), length(time2plot)*(FigProp.row-1)+aja) % use a new subplot
    imagesc(Clone_2Dmaps_2plot,[1 latticeDim^2]) % clone IDs are displayed on 2D space
    colormap(myColClones); % each clone ID is mapped to a different color
    
    %% Set axis properties:
    % force no ticks & square display:
    ax = gca; axis(ax,'off'); axis square;
    
    % Draw y-axis label with entry legend only for first-time-point (most left) panels:
    if aja == 1; ylabel(ax,FigProp.Leyend); set(get(ax,'YLabel'),'Visible','on'); end
    
    % Draw x-axis labels only for headers (top) panels:
    if setHeader == 1
        set(ax,'xaxisLocation','top')
        xlabel(ax, sprintf('%.1f',time2plot(aja))); set(get(ax,'XLabel'),'Visible','on')
    end
end

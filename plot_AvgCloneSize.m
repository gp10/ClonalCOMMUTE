function plot_AvgCloneSize(timepoints,avgCloneSize,showCI,FigProp)
%% Plots the average clone size over time
% The time course of the average clone size is plotted with or without
% confidence bounds (if provided).

% from Colom et al, 2020

%% Input:
% timepoints: vector of time points (expressed in weeks)
% avgCloneSize: structure containing average clone size (mandatory) and average clone size ± confidence bounds (optional)
    % struct{avg==..., ci95up==..., ci95dn==...}
        % avg: vector containing the average clone size at the different time points
        % ci95up: vector with upper confidence bound (same structure as avg)
        % ci95dn: vector with lower confidence bound (same structure as avg)
% showCI: calculate plausible intervals on model outcome (based on experimental sampling error) ( 0=NO | 1=YES )
% FigProp: structure containing output display properties
    % struct{Color=='r', Leyend=='CTL'}
        % Color: color used for plotting
        % Leyend: legend entry

%% Example:
% timepoints = [1:10];
% avgCloneSize.avg = 1:10; avgCloneSize.ci95up = [1:10]*1.1; avgCloneSize.ci95dn = [1:10]*0.9; 
% showCI = 1;
% FigProp.Color = 'r'; FigProp.Leyend = 'CTL';
% plot_AvgCloneSize(timepoints,avgCloneSize,showCI,FigProp);

%% Figure panel selection:
% Make handle to a new figure if this does not exist yet:
if isempty(findobj( 'Type', 'Figure', 'Name', 'avgCloneSizePlot'))
    % makes new figure:
    figure()
    set(gcf,'Name','avgCloneSizePlot');
else
    % retrieves existing one:
    figure(findobj( 'Type', 'Figure', 'Name', 'avgCloneSizePlot'))
end

%% Plot average clone size over time:
hold on
plot(timepoints.*7./365.*12,avgCloneSize.avg,'Color',FigProp.Color)
if showCI == 1
    plot(timepoints.*7./365.*12,avgCloneSize.ci95up,'Color',FigProp.Color)
    plot(timepoints.*7./365.*12,avgCloneSize.ci95dn,'Color',FigProp.Color)
end
xlabel('Time (months)'); ylabel('Avg. clone size')

%% Legend - retrieve preexisting legend and/or add new entry:
old_legend=findobj(gcf, 'Type', 'Legend');
if isempty(old_legend) % legend is empty; create one
    legend(FigProp.Leyend)
else
    legend([old_legend.String,FigProp.Leyend])
end

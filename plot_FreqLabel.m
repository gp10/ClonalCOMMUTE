function plot_FreqLabel(timepoints,FreqLabel,FigProp)
%% Plots the fraction of labelled basal cells over time
% The time course of the percentage of labelled basal cells is plotted.

% from Colom et al, 2020

%% Input:
% timepoints: vector of time points (expressed in weeks)
% FreqLabel: vector with the fraction of labelled basal cells at different time points 
% FigProp: structure containing output display properties
    % struct{Color=='r', Leyend=='CTL'}
        % Color: color used for plotting
        % Leyend: legend entry

%% Example:
% timepoints = [1:10];
% FreqLabel = ones(1,10); 
% FigProp.Color = 'r'; FigProp.Leyend = 'CTL';
% plot_FreqLabel(timepoints,FreqLabel,FigProp);

%% Figure panel selection:
% Make handle to a new figure if this does not exist yet:
if isempty(findobj( 'Type', 'Figure', 'Name', 'FreqLabel'))
    % makes new figure:
    figure()
    set(gcf,'Name','FreqLabel');
else
    % retrieves existing one:
    figure(findobj( 'Type', 'Figure', 'Name', 'FreqLabel'))
end

%% Plot fraction of labelled cells over time:
hold on
plot(timepoints.*7./365.*12, FreqLabel ,'Color',FigProp.Color)
ylim([0 10*max(FreqLabel)])
xlabel('Time (months)'); ylabel('% labelled cells')

%% Legend - retrieve preexisting legend and/or add new entry:
old_legend=findobj(gcf, 'Type', 'Legend');
if isempty(old_legend) % legend is empty; create one
    legend(FigProp.Leyend)
else
    legend([old_legend.String,FigProp.Leyend])
end

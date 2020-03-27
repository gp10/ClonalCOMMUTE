function plot_FreqMut(timepoints,FreqMut,FigProp)
%% Plots the fraction of mutant basal cells over time
% The time course of the percentage of mutant basal cells is plotted.

% from Colom et al, 2020

%% Input:
% timepoints: vector of time points (expressed in weeks)
% FreqMut: vector with the fraction of mutant basal cells at different time points 
% FigProp: structure containing output display properties
    % struct{Color=='r', Leyend=='CTL'}
        % Color: color used for plotting
        % Leyend: legend entry

%% Example:
% timepoints = [1:10];
% FreqMut = ones(1,10); 
% FigProp.Color = 'r'; FigProp.Leyend = 'CTL';
% plot_FreqMut(timepoints,FreqMut,FigProp);

%% Figure panel selection:
% Make handle to a new figure if this does not exist yet:
if isempty(findobj( 'Type', 'Figure', 'Name', 'FreqMut'))
    % makes new figure:
    figure()
    set(gcf,'Name','FreqMut');
else
    % retrieves existing one:
    figure(findobj( 'Type', 'Figure', 'Name', 'FreqMut'))
end

%% Plot fraction of mutant cells over time:
hold on
plot(timepoints.*7./365.*12, FreqMut ,'Color',FigProp.Color)
xlabel('Time (months)'); ylabel('% mutant cells')

%% Legend - retrieve preexisting legend and/or add new entry:
old_legend=findobj(gcf, 'Type', 'Legend');
if isempty(old_legend) % legend is empty; create one
    legend(FigProp.Leyend)
else
    legend([old_legend.String,FigProp.Leyend])
end

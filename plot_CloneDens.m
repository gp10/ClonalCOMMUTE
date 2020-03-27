function plot_CloneDens(timepoints,NClones,showCI,FigProp)
%% Plots the clone density over time
% The time course of the number surviving clones (normalized to initial
% time point) is plotted with or without confidence bounds (if provided).

% from Colom et al, 2020

%% Input:
% timepoints: vector of time points (expressed in weeks)
% CloneDens: structure containing clone density (mandatory) and clone density ± confidence bounds (optional)
    % struct{avg==..., ci95up==..., ci95dn==...}
        % avg: vector containing the average clone density at the different time points
        % ci95up: vector with upper confidence bound (same structure as avg)
        % ci95dn: vector with lower confidence bound (same structure as avg)
% showCI: calculate plausible intervals on model outcome (based on experimental sampling error) ( 0=NO | 1=YES )
% FigProp: structure containing output display properties
    % struct{Color=='r', Leyend=='CTL'}
        % Color: color used for plotting
        % Leyend: legend entry

%% Example:
% timepoints = [1:10];
% NClones.avg = exp(-1*[1:10]); NClones.ci95up = exp(-1*[1:10])*1.1; NClones.ci95dn = exp(-1*[1:10])*0.9; 
% showCI = 1;
% FigProp.Color = 'r'; FigProp.Leyend = 'CTL';
% plot_CloneDens(timepoints,NClones,showCI,FigProp);

%% Figure panel selection:
% Make handle to a new figure if this does not exist yet:
if isempty(findobj( 'Type', 'Figure', 'Name', 'CloneDensPlot'))
    % makes new figure:
    figure()
    set(gcf,'Name','CloneDensPlot');
else
    % retrieves existing one:
    figure(findobj( 'Type', 'Figure', 'Name', 'CloneDensPlot'))
end

%% Plot clone density over time:
hold on
plot(timepoints.*7./365.*12, NClones.avg ./ NClones.avg(1).*100 ,'Color',FigProp.Color)
if showCI == 1
    plot(timepoints.*7./365.*12, NClones.ci95up ./ NClones.ci95up(1).*100 ,'Color',FigProp.Color)
    plot(timepoints.*7./365.*12, NClones.ci95dn ./ NClones.ci95dn(1).*100 ,'Color',FigProp.Color)
end
set(gca,'YScale','log')
xlabel('Time (months)'); ylabel('Clone density (%) (rel. to initial)')

%% Legend - retrieve preexisting legend and/or add new entry:
old_legend=findobj(gcf, 'Type', 'Legend');
if isempty(old_legend) % legend is empty; create one
    legend(FigProp.Leyend)
else
    legend([old_legend.String,FigProp.Leyend])
end

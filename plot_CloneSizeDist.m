function plot_CloneSizeDist(timepoints,time2plot,CUM_nfreq,showCI,avgCloneSize,FigProp)
%% Plots normalized, cumulative clone size distributions over time
% It plots the cumulative distribution of clone sizes at time points of
% interest. Clone sizes (n) are normalized by the average clone size (<n>)
% at each specific time point.

% from Colom et al, 2020

%% Input:
% timepoints: vector of time points when clone size frequencies have been calculated (expressed in weeks)
% time2plot: vector of time points of interest to plot data from (expressed in weeks)
% CUM_nfreq: structure containing clone size frequencies (mandatory) and clone size frequencies ± confidence bounds (optional)
    % struct{full==..., centre==..., ci95up==..., ci95dn==...}
        % full: cell array of normalized cumulative clone size frequencies at given time points
        % centre: mid of 95 bound on normalized cumulative clone size frequencies (sampling error)
        % ci95up: upper 95% bound on normalized cumulative clone size frequencies (sampling error)
        % ci95dn: bottom 95% bound on normalized cumulative clone size frequencies (sampling error)
% showCI: calculate plausible intervals on model outcome (based on experimental sampling error) ( 0=NO | 1=YES )
% avgCloneSize: vector of average clone sizes at different time points
% FigProp: structure containing output display properties
    % struct{Color=='r', Leyend=='CTL'}
        % Color: color used for plotting
        % Leyend: legend entry

%% Example:
% timepoints = [1:3];
% time2plot = [2,3];
% CUM_nfreq.full{1,1} = fliplr(log10([1:10]))'; CUM_nfreq.full{1,2} = fliplr(log10([1:0.1:10]))'; CUM_nfreq.full{1,3} = fliplr(log10([1:0.01:10]))';
% CUM_nfreq.centre = {}; CUM_nfreq.ci95up = {}; CUM_nfreq.ci95dn = {};
% showCI = 0;
% avgCloneSize = [10,50,120];
% FigProp.Color = 'r'; FigProp.Leyend = 'CTL';
% plot_CloneSizeDist(timepoints,time2plot,CUM_nfreq,showCI,avgCloneSize,FigProp);

%% Figure panel selection:
% Make handle to a new figure if this does not exist yet:
if isempty(findobj( 'Type', 'Figure', 'Name', 'CloneSizeDist'))
    % makes new figure (with white background color):
    figure()
    set(gcf,'Name','CloneSizeDist');
else
    % retrieves existing one:
    figure(findobj( 'Type', 'Figure', 'Name', 'CloneSizeDist'))
end

%% Plot normalized, cumulative clone size distributions at the desired time points:
for aja = 1:length(time2plot)
    myTime = find(timepoints >= time2plot(aja),1); % position in the time vector
    subplot(length(time2plot),length(time2plot),aja) % generate time-specific subplots
    hold on
    % cumulative clone size distributions are normalized by the corresponding average clone size at each time point:
    plot([0:size(CUM_nfreq.full{1,myTime},1)-1]./avgCloneSize(myTime), CUM_nfreq.full{1,myTime},...
        'Color',FigProp.Color,'Marker','none','LineStyle','-','LineWidth',2)
    if showCI == 1
        plot([0:size(CUM_nfreq.centre{1,myTime},1)-1]./avgCloneSize(myTime), CUM_nfreq.ci95up{1,myTime},'Color',FigProp.Color)
        plot([0:size(CUM_nfreq.centre{1,myTime},1)-1]./avgCloneSize(myTime), CUM_nfreq.ci95dn{1,myTime},'Color',FigProp.Color)
    end
    
    % Settings:
    set(gca,'YScale','log')
    ylim([1E-3 1]); set(gca,'YTick',[0.001 0.01 0.1 1])
    xlim([0 10]); set(gca,'XTick',[0:5:10])
    box on; grid on;
    %plot([0:1:20],exp(-[0:1:20]),'k-')
    xlabel('n / <n>'); ylabel('Frequency')
end

%% Legend - retrieve preexisting legend and/or add new entry:
old_legend=findobj(gcf, 'Type', 'Legend');
if isempty(old_legend) % legend is empty; create one
    legend(FigProp.Leyend)
else
    legend([old_legend.String,FigProp.Leyend])
end

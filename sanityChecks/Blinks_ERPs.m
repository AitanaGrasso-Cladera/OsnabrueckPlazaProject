%% This code generates ERP plots using blink onsets (eyes starting to close)

clc; 
close all; 

%%     
% load eeglab and xdf plugin
addpath('');
addpath('');
eeglab; close; 

% get files
dataFolder = '';

% list participant folders
tmp = dir(fullfile(dataFolder, 'P*'));

participants = struct();
inx = 1;

for k = 1:numel(tmp)

    participants(inx).folder = fullfile(tmp(k).folder, tmp(k).name);
    participants(inx).name     = tmp(k).name;   
    participants(inx).id     = tmp(k).name(2:end);
    inx = inx + 1;
end

allERP = [];   % will store Cz ERP for each recording for grand erp plot

% Load xdf, triggers, hardware triggers
for sub = 1:numel(participants)

    participantFolder = participants(sub).folder; 
  
    % Load epoched data
    EEG = pop_loadset(sprintf('5a_epochedBlinkRejection_%s.set',participants(sub).id),fullfile(participantFolder));
    
    % ERP calculation
    erp = mean(EEG.data, 3);   % channels × time
    Cz = find(strcmpi({EEG.chanlocs.labels}, 'Cz'));
    % Cz ERP
    czERP = erp(Cz, :);     % for Grand average erp 

    allERP = [allERP; czERP];   % recordings × time    
    
    %% ERP images with sorted blinks

    nEvents = length(EEG.event);
    blinkDurations = zeros(nEvents, 1);
    saccadeDurations = zeros(nEvents, 1);
    for i = 1:nEvents
       dur_str = EEG.event(i).givenDurations; % char 
       sacc_str = EEG.event(i).saccadeDurations; 
       EEG.event(i).givenDurations = str2double(dur_str); % convert to numeric 
       EEG.event(i).saccadeDurations = str2double(sacc_str); 
    end

    EEG = eeg_checkset(EEG, 'eventconsistency');
    
    figure('Position', [500 1200 1200 500]);  
     
    % sorted by saccade duration preceeding blinks
    subplot(1,2,1);
    erpimage( mean(EEG.data([11], :),1), eeg_getepochevent(EEG, {'empty'}, [], 'saccadeDurations'), ...
        linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts), '', 10, 1 ,'yerplabel','\muV','erp','off','colorbar','on',...
        'replace_ties', 'off');
    title('Sorted by closest saccade durations before the blink');

    % sorted by blink duration
    subplot(1,2,2);
    erpimage( mean(EEG.data([11], :),1), eeg_getepochevent(EEG, {'empty'}, [], 'givenDurations'), ...
        linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts), '', 10, 1 ,'yerplabel','\muV','erp','off','colorbar','on',...
        'replace_ties', 'off');
    title('Sorted  by blink durations');

    %exportgraphics(gcf, sprintf('blinks_sorted_%s.png', participants(sub).id), ...
      %'Resolution', 300, 'Width', 12, 'Height', 4, 'Units', 'inches')


    %% ERP images with trials sorted by time (default)
    figure;
    erpimage( mean(EEG.data([11], :),1), ones(1, EEG.trials)*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts), 'Blink ERP, sorted by trial onset', 10, 1 ,'yerplabel','\muV','erp','on','colorbar','on',...
        'replace_ties', 'off');

    %exportgraphics(gcf, sprintf('blinks_%s.png', participants(sub).id), ...
      %'Resolution', 300, 'Width', 8, 'Height', 6, 'Units', 'inches')
    %}
end

%% Grand average ERP for Cz with CI of 95%
grandERP = mean(allERP, 1);   % average across recordings
times = EEG.times;     

% number of recordings
nRec = size(allERP, 1);

% standard deviation
std_allERP = std(allERP, 0, 1);

% 95% CI using t-distribution (2 tails
% (makes 0.975 the upper limit)
tval = tinv(0.975, nRec-1);
CI = tval * (std_allERP./ sqrt(nRec));

figure; hold on
fill([times fliplr(times)], ...
     [grandERP+CI fliplr(grandERP-CI)], ...
     [0.8 0.8 1], 'EdgeColor','none');
plot(times, grandERP, 'b', 'LineWidth', 2);

xline(0,'--k');
yline(0, '--k');
xlabel('Time (ms)');
ylabel('Amplitude (\muV)');
title('Grand ERP at Cz (95% CI)');
grid on

%% grand ERP with all ERPs per recording

figure; hold on
colors = jet(nRec) * 0.90;  % colour brightness
legendHandles = gobjects(nRec, 1);  %  legend handles
legendLabels = cell(nRec, 1);       % labels

for erp = 1:nRec
    h = plot(times, allERP(erp, :), ...
             'Color', colors(erp, :), ...
             'LineWidth', 1.5);  % lines width
    legendHandles(erp) = h;  
    % store participant id
    legendLabels{erp} = participants(erp).id; 
end

plot(times, grandERP, 'k', 'LineWidth', 3);
legend(legendHandles, legendLabels, 'Location', 'best');

xline(0, '--k');
yline(0, '--k');
xlabel('Time (ms)');
ylabel('Amplitude (\muV)');
title('Individual Cz ERPs with Grand Average');


grid on
hold off;

%%
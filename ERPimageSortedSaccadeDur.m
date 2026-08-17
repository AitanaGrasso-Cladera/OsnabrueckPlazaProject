%% Script to generate plots sorted by saccade duration and plots to explore data
%% History
% 25.02.2025 = AGC: The scripts sorts fixation ERPs based on the duration
% of the preceeding saccade. Then, it plots the data for a specific channel
% of interest (for all subjects).
%%
projectFolder = [];
dataFolder = [projectFolder,filesep,'Data'];
addpath('[]/eeglab2026.0.0');
workingFolder = [];
figureFolder = [projectFolder,filesep,'figures'];
%%
tmp = dir(fullfile(dataFolder));
participants = [];
inx = 1;

for pId = 1:size(tmp,1)
    if tmp(pId).name(1) == '.'|| contains(tmp(pId).name,'xdf') || ~contains(tmp(pId).name,'preproc')
        continue
    else
        participants(inx).name = tmp(pId).name(9:end);
        participants(inx).folder = tmp(pId).folder;
        participants(inx).date = tmp(pId).date;
        inx = inx+1;
    end
end

clear tmp
%%
eeglab
close
%%
for sub = 1:length(participants)

    participantFolder = [dataFolder,filesep,'preproc_',participants(sub).name];

    EEGfix = pop_loadset(sprintf('4a_epochedFixationPlot_%s.set',participants(sub).name),fullfile(participantFolder));

    numEpochs = length(EEGfix.epoch);
    latency_array = NaN(numEpochs,1);  % Preallocate for speed
    type_array = cell(numEpochs,1);
    saccade_dur = cell(numEpochs,1);
    saccade_onset = cell(numEpochs,1);

    for i = 1:numEpochs
        % Extract latency cell and convert to array
        latencies = cell2mat(EEGfix.epoch(i).eventlatency);

        % Find index where latency is 0
        zero_idx = find(latencies == 0, 1);

        if ~isempty(zero_idx)
            % Extract corresponding values from other fields
            latency_array(i) = latencies(zero_idx);
            type_array{i} = EEGfix.epoch(i).eventtype{zero_idx};
            saccade_dur{i} =  -str2double(EEGfix.epoch(i).eventgivenDurations{zero_idx});
        end
    end

    event = EEGfix.event(:,1:numEpochs);
    [event.saccade_dur] = deal(saccade_dur{:});

    % Sort saccade durations in ascending order
    chn = 11;
    data = squeeze(EEGfix.data(chn,:,:));
    times = EEGfix.times;
    event = EEGfix.event(:,1:numEpochs);
    [event.saccade_dur] = deal(saccade_dur{:});
    % Sort saccade durations in ascendent order
    saccadeDurations = [event.saccade_dur];
    [sortedSaccades,sortedIdxs] = sort(saccadeDurations,'ascend');
    % Sort data based on saccade duration
    sortedData = data(:,sortedIdxs);
    % Parameters
    binSize = 30;
    nTrials = size(sortedData, 2);
    % assuming trials are columns
    nBins = floor(nTrials / binSize);
    % Preallocate
    binnedData = zeros(size(sortedData, 1), nBins);
    % Loop over bins
    saccadeOnsets = zeros(1, nBins);
    % preallocate
    for b = 1:nBins
        idxStart = (b - 1) * binSize + 1;
        idxEnd = b * binSize;
        binnedData(:, b) = mean(sortedData(:, idxStart:idxEnd), 2);
        saccadeOnsets(b) = mean(sortedSaccades(idxStart:idxEnd));
        % scalar per bin
    end

    f = figure;
    f.Position = [575,45,1043,860];
    subplot(4,1,1)
    plot(times, squeeze(mean(data,2)), 'LineWidth', 1.5,'Color','b')
    xlim([min(times), max(times)])
    xline(0, 'k-', 'LineWidth', 2);
    xlabel('Time (ms)','FontSize',14);
    ylabel('Amplitude (mV)','FontSize',14)
    ylim([min(squeeze(mean(data,2)))-0.5 max(squeeze(mean(data,2)))+0.5])
    title('Fixation Onset','FontSize',14)
    set(gca, 'FontSize', 14) % axis ticks
    subplot(4,1,2)
    plot(times, squeeze(mean(data2,2)), 'LineWidth', 1.5,'Color','r')
    xlim([min(times), max(times)])
    xline(0, 'k-', 'LineWidth', 2);
    xlabel('Time (ms)','FontSize',14);
    ylim([min(squeeze(mean(data2,2)))-0.5 max(squeeze(mean(data2,2)))+0.5])
    ylabel('Amplitude (mV)','FontSize',14)
    title('Saccade Onset','FontSize',14)
    set(gca, 'FontSize', 14) % axis ticks
    subplot(4,1,3)
    plot(times, squeeze(mean(data3,2)), 'LineWidth', 1.5,'Color','k')
    xlim([min(times), max(times)])
    xline(0, 'k-', 'LineWidth', 2);
    xlabel('Time (ms)','FontSize',14);
    ylim([min(squeeze(mean(data3,2)))-0.5 max(squeeze(mean(data3,2)))+0.5])
    ylabel('Amplitude (mV)','FontSize',14)
    title('Blink Onset','FontSize',14)
    set(gca, 'FontSize', 14) % axis ticks
    subplot(4,1,4)
    plot(times, squeeze(mean(data4,2)), 'LineWidth', 1.5,'Color','k')
    xlim([min(times), max(times)])
    xline(0, 'k-', 'LineWidth', 2);
    xlabel('Time (ms)','FontSize',14);
    ylim([min(squeeze(mean(data3,2)))-0.5 max(squeeze(mean(data3,2)))+0.5])
    ylabel('Amplitude (mV)','FontSize',14)
    title('Blink Offset','FontSize',14)
    set(gca, 'FontSize', 14) % axis ticks
    fileName = ['ERPwaveComparison_',participants(sub).name,'.png'];
    %exportgraphics(gcf, fullfile(participantFolder,fileName), 'Resolution', 600);
    %close

    ALLDATAsacc{sub} = EEGsacc.data;
    ALLDATAfix{sub} = EEGfix.data;
    ALLDATAblink{sub} = EEGblink.data;

    % Save participant data
    binnedDataALL{sub} = binnedData;
    saccadeOnsetsALL{sub} = saccadeOnsets;
end
save(fullfile(dataFolder,sprintf('binnedDataSacc.mat')),'binnedDataSaccALL','saccadeOnsetsALL');



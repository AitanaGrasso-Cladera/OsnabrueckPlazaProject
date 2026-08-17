%% Script to generate plots sorted by blink duration and plots to explore data
%% History
% April.2026 = AGC: The scripts sorts blink onset ERPs based on the duration
% of the blink event. Then, it plots the data for a specific channel
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

    EEGblink = pop_loadset(sprintf('4a_epochedBlinkOffset_%s.set',participants(sub).name),fullfile(participantFolder));

    numEpochs = length(EEGblink.epoch);
    latency_array = NaN(numEpochs,1);  % Preallocate for speed
    type_array = cell(numEpochs,1);
    blink_dur = cell(numEpochs,1);
    saccade_onset = cell(numEpochs,1);

    for i = 1:numEpochs
        % Extract latency cell and convert to array
        latencies = cell2mat(EEGblink.epoch(i).eventlatency);

        % Find index where latency is 0
        zero_idx = find(latencies == 0, 1);

        if ~isempty(zero_idx)
            % Extract corresponding values from other fields
            latency_array(i) = latencies(zero_idx);
            type_array{i} = EEGblink.epoch(i).eventtype{zero_idx};
            blink_dur{i} =  -str2double(EEGblink.epoch(i).eventgivenDurations{zero_idx});
        end
    end

    event = EEGblink.event(:,1:numEpochs);
    [event.blink_dur] = deal(blink_dur{:});

    % Sort blink durations in ascending order
    chn = 11;
    data = squeeze(EEGblink.data(chn,:,:));
    times = EEGblink.times;
    event = EEGblink.event(:,1:numEpochs);
    [event.blink_dur] = deal(blink_dur{:});
    % Sort blink durations in ascendent order
    blinkDurations = [event.blink_dur];
    [sortedBlinks,sortedIdxs] = sort(blinkDurations,'ascend');
    % Sort data based on blink duration
    sortedData = data(:,sortedIdxs);
    % Parameters
    binSize = 30;
    nTrials = size(sortedData, 2);
    % assuming trials are columns
    nBins = floor(nTrials / binSize);
    % Preallocate
    binnedData = zeros(size(sortedData, 1), nBins);
    % Loop over bins
    blinkOnsets = zeros(1, nBins);
    % preallocate
    for b = 1:nBins
        idxStart = (b - 1) * binSize + 1;
        idxEnd = b * binSize;
        binnedData(:, b) = mean(sortedData(:, idxStart:idxEnd), 2);
        blinkOnsets(b) = mean(sortedSaccades(idxStart:idxEnd));
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

    % Save participant data
    binnedDataBlinkALL{sub} = binnedData;
    blinkOnsetsALL{sub} = blinkOnsets;
end
save(fullfile(dataFolder,sprintf('binnedDataBlink.mat')),'binnedDataBlinkALL','blinkOnsetsALL');



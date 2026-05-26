%% Unfold for Event Related Potentials
% This script was developed by Aitana Grasso-Cladera and Debora Nolte
%% History
% 03.06.25 - DN - AGC = Created the first version to run Unfold based in
% Nolte et al., 2025 (eNeuro) and Grasso-Cladera et al., 2025 (Brain
% Sciences)
% 07.06.25 - AGC = Following Ladouce et al., 2023, we removed aberrant
% activity (threshold of 5 std around the median voltage activity recorded
% across all epochs) were discarded.
% --.12.25 - AGC = Code now applies a 2.5 threshold calculation based on
% IQR to remove bad trials.
% 02.02.26 - AGC = Epoching and removing trials for blink onset ERPs.
% 12.02.26 - AGC = Added the IQR rejection procedure for cleaning the data.
%% Path definition
projectFolder = [];
dataFolder = [projectFolder,filesep,'Data'];
addpath(genpath('[]/eeglab2026.0.0'));
workingFolder = [];
addpath('[]/unfold')
init_unfold
addpath('[]/unfold/gramm')
%%
eeglab;
close
%%
tmp = dir(fullfile(dataFolder));
participants = [];
inx = 1;

for pId = 1:size(tmp,1)
    if tmp(pId).name(1) == '.'|| ~contains(tmp(pId).name,'preproc_')
        continue
    else
        participants(inx).name = tmp(pId).name;
        participants(inx).folder = tmp(pId).folder;
        participants(inx).date = tmp(pId).date;
        inx = inx+1;
    end
end

clear tmp
%%
cnt = 1;
for sub = 1:length(participants)
  
    participantFolder = [dataFolder,filesep,participants(sub).name];
    % add this new folder to the savedata path, so where the intermediate steps
    % will be saved
    savedata = [dataFolder,filesep,'unfold'];
    if ~exist(savedata,'dir')
        mkdir(savedata)
    end
    %% Load the data
    EEG = pop_loadset(sprintf('2a_cleanDataChannels_noRejectionAll_%s.set',participants(sub).name(9:end)),fullfile([dataFolder,filesep,participants(sub).name]));
    %% Interpolating missing channels
    % Get all channels that need to be interpolated
    EEGchan = pop_loadset(sprintf('1a_triggersFilteredAll_%s.set',participants(sub).name(9:end)),fullfile([dataFolder,filesep,participants(sub).name]));

    fullChanlocs = EEGchan.chanlocs; % Used for data cleaning and interpolation
    clear EEGchan
    EEG = pop_interp(EEG,fullChanlocs,'spherical');
    EEG = eeg_checkset(EEG);
    %% Prepare the event file
    % Remove boudary events (for now)
    boundaryEvents = false(size(EEG.event));
    for j = 1:length(EEG.event)
        if contains(EEG.event(j).type, 'boundary')
            boundaryEvents(j) = true;
        end
    end

    EEG.event(boundaryEvents) = [];

    % Remove invalid fixations and saccades
    invalidEvents = false(size(EEG.event));
    for i = 1:length(EEG.event)
        if contains(EEG.event(i).status, 'invalid')
            invalidEvents(i) = true;
        end
    end
    EEG.event(invalidEvents) = [];

    EEG1 = EEG;
    %% Create design matrix
    cfgDesign               = [];
    % all events we have are saccades, blinks and fixations:
    cfgDesign.eventtypes    = {'saccadeOnset','blinkOnset','blinkOffset','fixationOnset'};
    cfgDesign.codingschema  = 'effects';
    cfgDesign.formula       = {'y ~ 1','y ~ 1','y ~ 1','y ~ 1'};
    
    EEG                    = uf_designmat(EEG,cfgDesign);
    %% Timeshift
    cfgTimeshift            = [];
    cfgTimeshift.timelimits = [-.5 .7];
    EEG                    = uf_timeexpandDesignmat(EEG,cfgTimeshift);
    %% Reject noisy segments
    % % we only do this, if we have noisy segments to reject
    if isfile(fullfile([dataFolder,filesep,participants(sub).name],sprintf('removed_intervalsCleanRawDataNEW_%s.mat',participants(sub).name(9:end))))
        load(fullfile([dataFolder,filesep,participants(sub).name],sprintf('removed_intervalsCleanRawDataNEW_%s.mat',participants(sub).name(9:end))));
        % Check if we are removing the blink completely (on + off), if not,
        % do it
        eventTypes = {EEG.event.type};
        onsetIdx  = find(strcmp(eventTypes,'blinkOnset'));
        offsetIdx = find(strcmp(eventTypes,'blinkOffset'));

        newWinrej = [];
        for i = 1:length(onsetIdx)
            onset_i   = onsetIdx(i);
            onset_lat = EEG.event(onset_i).latency;
            % Find first offset occurring after this onset
            offset_i = offsetIdx(find(offsetIdx > onset_i,1,'first'));
            if isempty(offset_i)
                continue
            end
            offset_lat = EEG.event(offset_i).latency;
            % Check whether onset or offset is already rejected
            onsetRejected  = any(onset_lat  >= tmprej(:,1) & onset_lat  <= tmprej(:,2));
            offsetRejected = any(offset_lat >= tmprej(:,1) & offset_lat <= tmprej(:,2));
            % If either one is rejected → reject the full blink interval
            if onsetRejected || offsetRejected
                newWinrej = [newWinrej; onset_lat offset_lat];
            end
        end
        % Merge with the matrix coming from clean raw data, sort them and
        % merge overlapping windows
        tmprej = [tmprej; newWinrej];
        if ~isempty(tmprej)
            tmprej = sortrows(tmprej);
            merged = tmprej(1,:);
            for i = 2:size(tmprej,1)
                if tmprej(i,1) <= merged(end,2)
                    merged(end,2) = max(merged(end,2), tmprej(i,2));
                else
                    merged = [merged; tmprej(i,:)];
                end
            end
            tmprej = merged;
        end

        EEG                = uf_continuousArtifactExclude(EEG,'winrej',tmprej);
    end
    %% Detect noisy segments after clean raw data removal
    fs = 500;
    winLength = round(1*fs); % 500 ms window
    stepSize = round(0.5*fs); % 0.05 is 250 ms step if 1 = fully slidding; if windowLenght no overlapping
    % 0.5 is half the window
    nSamples = size(EEG.data,2);
    nChans   = size(EEG.data,1);

    winStarts = 1:stepSize:(nSamples - winLength + 1);
    nWins = length(winStarts);
    % IQR per Window - similar to trial level before
    IQRwindow = nan(1,nWins);

    for i = 1:nWins
        windowIDX = winStarts(i):(winStarts(i)+winLength-1);
        windowData = EEG.data(:,windowIDX);
        IQRwindow(i) = iqr(windowData(:));
    end

    % Global threshold
    medianIQRwins = median(IQRwindow);
    generalIQRwins = iqr(IQRwindow);

    lowerThreshold = medianIQRwins - (generalIQRwins*3.5);
    upperThreshold = medianIQRwins + (generalIQRwins*3.5);

    % Channel level IQR per window
    IQRchannels = nan(nChans, nWins);
    bad_windows = false(1,nWins);

    for i = 1:nWins
        idx = winStarts(i):(winStarts(i)+winLength-1);
        windowData = EEG.data(:,idx);
        for ch = 1:nChans
            IQRchannels(ch,i) = iqr(windowData(ch,:));
        end
        % Mark window as bad if ANY channel exceeds threshold
        if any(IQRchannels(:,i) > upperThreshold) || ...
                any(IQRchannels(:,i) < lowerThreshold)
            bad_windows(i) = true;
        end
    end

    badIdx = find(bad_windows); % indices of bad windows
    if ~isempty(badIdx)
        % Convert window indices to sample indices
        winStartSamples = winStarts(badIdx);
        winEndSamples   = winStartSamples + winLength - 1;
        % Merge overlapping windows (important because windows overlap!)
        starts = winStartSamples(1);
        ends   = [];
        currentEnd = winEndSamples(1);
        for i = 2:length(winStartSamples)
            if winStartSamples(i) <= currentEnd + 1
                % Overlapping or touching → extend current segment
                currentEnd = max(currentEnd, winEndSamples(i));
            else
                % New segment
                ends = [ends, currentEnd];
                starts = [starts, winStartSamples(i)];
                currentEnd = winEndSamples(i);
            end
        end
        % Add last segment
        ends = [ends, currentEnd];
        bad_windowsTimes = [starts; ends]';
    else
        bad_windowsTimes = [];
    end

    % Check if we are removing the blink completely (on + off), if not,
        % do it
        eventTypes = {EEG.event.type};
        onsetIdx  = find(strcmp(eventTypes,'blinkOnset'));
        offsetIdx = find(strcmp(eventTypes,'blinkOffset'));

        newWinrej = [];
        for i = 1:length(onsetIdx)
            onset_i   = onsetIdx(i);
            onset_lat = EEG.event(onset_i).latency;
            % Find first offset occurring after this onset
            offset_i = offsetIdx(find(offsetIdx > onset_i,1,'first'));
            if isempty(offset_i)
                continue
            end
            offset_lat = EEG.event(offset_i).latency;
            % Check whether onset or offset is already rejected
            onsetRejected  = any(onset_lat  >= bad_windowsTimes(:,1) & onset_lat  <= bad_windowsTimes(:,2));
            offsetRejected = any(offset_lat >= bad_windowsTimes(:,1) & offset_lat <= bad_windowsTimes(:,2));
            % If either one is rejected → reject the full blink interval
            if onsetRejected || offsetRejected
                newWinrej = [newWinrej; onset_lat offset_lat];
            end
        end
        % Merge with the matrix coming from clean raw data, sort them and
        % merge overlapping windows
        bad_windowsTimes = [bad_windowsTimes; newWinrej];
        merged = [];
        if ~isempty(bad_windowsTimes)
            bad_windowsTimes = sortrows(bad_windowsTimes);
            merged = bad_windowsTimes(1,:);
            for i = 2:size(bad_windowsTimes,1)
                if bad_windowsTimes(i,1) <= merged(end,2)
                    merged(end,2) = max(merged(end,2), bad_windowsTimes(i,2));
                else
                    merged = [merged; bad_windowsTimes(i,:)];
                end
            end
            bad_windowsTimes = merged;
        end

    save(fullfile(participantFolder,sprintf('removed_IQRcontinuousFinalMatchingBlinks_%s.mat',participants(sub).name(9:end))),'bad_windows','bad_windowsTimes');
    %% Reject noisy segments
    % % we only do this, if we have noisy segments to reject
    if ~isempty(bad_windowsTimes)
        EEG                = uf_continuousArtifactExclude(EEG,'winrej',bad_windowsTimes);
    end
    %% Fit the model
    % fit the model:
    % the parameters are taken from Gert et al., 2022
    cfgFit                  = [];
    cfgFit.precondition     = 1;
    cfgFit.lsmriterations   = 1500; % steps iterative solver should reach
    cfgFit.channel          = 1:length(EEG.chanlocs); % all channels
    EEG                    = uf_glmfit(EEG, cfgFit);
    %% Make a massive uni-variate fit without de-convolution (Gert et al., 2022)
    EEGepoch = uf_epoch(EEG,'timelimits',cfgTimeshift.timelimits);

    EEGepoch = uf_glmfit_nodc(EEGepoch);
    %% Get the betas
    % results condensed in new structure
    ufresult                = uf_condense(EEG);
    ufresultEp              = uf_condense(EEGepoch);

    ufresultEp = uf_predictContinuous(ufresultEp); 
    ufresultEp = uf_addmarginal(ufresultEp);

    paramNames={ufresultEp.param.name};
    [~,paramSaccade]=find(ismember(paramNames,'(Intercept)'));
    [~,paramBlinkOnset]=find(ismember(paramNames,'2_(Intercept)'));
    [~,paramBlinkOnffset]=find(ismember(paramNames,'3_(Intercept)'));
    [~,paramFixation]=find(ismember(paramNames,'4_(Intercept)'));

    sacc(cnt,:,:) = ufresultEp.beta(:,:,paramSaccade);
    blinkOn(cnt,:,:)  = ufresultEp.beta(:,:,paramBlinkOnset);
    blinkOff(cnt,:,:)  = ufresultEp.beta(:,:,paramBlinkOnffset);
    fix(cnt,:,:)  = ufresultEp.beta(:,:,paramFixation);

    time = ufresultEp.times;

    cnt = cnt + 1; % to avoid having nans in the data

    %% Save the data
    save(fullfile(savedata,sprintf('betasUnfoldFinalMatchingBlinks_%s.mat',participants(sub).name)),'ufresult','ufresultEp');
end

save(fullfile(savedata,sprintf('betasUnfoldALL_overlapFinalMatchingBlinks.mat')),'sacc','blinkOn','fix','time','blinkOff');

%%
rmpath('[]/eeglab2026.0.0');
rmpath('[]/unfold')